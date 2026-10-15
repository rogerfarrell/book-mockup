/* book-mockup.c
 *
 * Copyright 2026 Roger
 *
 * SPDX-License-Identifier: GPL-3.0-or-later
 */

/**
 * book-mockup:
 *
 * Renders an upright book mockup at any rotation from flat front cover,
 * spine and back cover artwork, using a single GEGL graph.
 *
 * The code is organized in four strata. Each stratum depends only on the
 * strata above it, so any prefix of the file can be copied into another
 * project unchanged:
 *
 * 1. Math      — points, quads, homographies, convex hulls. No GEGL calls.
 * 2. Scene     — a pinhole camera, textured planar faces, a rotating
 *                six-sided box, lighting, and the floor shadow. No GEGL
 *                calls.
 * 3. Render    — generic GEGL graph helpers: load, shade, warp, mask,
 *                stack, save.
 * 4. App       — the book: options, command line, and main().
 *
 * Strata 1–3 know nothing about books. A box-shaped product (a software
 * box, a DVD case, a carton) needs only a new stratum 4.
 *
 * To reuse a stratum under another prefix, rename `bm_`, `Bm` and `BM_`
 * consistently, for example with:
 *
 * |[
 *   sed -e 's/\bbm_/xx_/g; s/\bBm/Xx/g; s/\bBM_/XX_/g'
 * ]|
 *
 * Coordinates: x grows to the right, y grows down (screen convention), and
 * z grows into the scene. y = 0 is eye level.
 *
 * Colors and shading are in linear light, as GEGL works internally.
 *
 * Build (macOS / Homebrew):
 *
 * |[
 *   brew install gegl pkgconf
 *   cc -std=c11 -O2 -Wall -o book-mockup book-mockup.c \
 *      $(pkg-config --cflags --libs gegl-0.4) -lm
 * ]|
 */

#include <math.h>
#include <stdlib.h>
#include <string.h>

#include <gegl.h>

#define BM_DEG_TO_RAD(d) ((d) * G_PI / 180.0)


/* Stratum 1: Math
 *
 * Plane geometry only. Depends on GLib types, <math.h> and the
 * `GeglMatrix3` struct (a plain 3×3 array, row-major: coeff[row][column]).
 */

#define BM_MAX_POLYGON_POINTS 16

typedef struct _BmPoint   BmPoint;
typedef struct _BmVec3    BmVec3;
typedef struct _BmQuad    BmQuad;
typedef struct _BmPolygon BmPolygon;

struct _BmPoint
{
  double x;
  double y;
};

struct _BmVec3
{
  double x;
  double y;
  double z;
};

/* Corners in order: top-left, top-right, bottom-right, bottom-left */
struct _BmQuad
{
  BmPoint p[4];
};

struct _BmPolygon
{
  BmPoint points[BM_MAX_POLYGON_POINTS];
  int     n_points;
};

/**
 * bm_quad_from_rect:
 * @width: rectangle width
 * @height: rectangle height
 *
 * Builds the quad of an axis-aligned rectangle at the origin, for use as
 * the source side of a homography.
 *
 * Returns: the quad
 */
static BmQuad
bm_quad_from_rect (double width,
                   double height)
{
  return (BmQuad) { .p = { { 0, 0 }, { width, 0 },
                           { width, height }, { 0, height } } };
}

/**
 * bm_homography_from_quads:
 * @source: the source quad
 * @dest: the destination quad
 * @matrix: (out caller-allocates): return location for the matrix
 *
 * Solves the projective transform that maps @source onto @dest
 * (four-point DLT, 8×8 Gaussian elimination with partial pivoting).
 *
 * Returns: %FALSE if the quads are degenerate
 */
static gboolean
bm_homography_from_quads (const BmQuad *source,
                          const BmQuad *dest,
                          GeglMatrix3  *matrix)
{
  double a[8][9];
  double h[8];
  int i, row, col;

  /* For each correspondence (x, y) → (u, v), with h33 = 1:
   *   h11 x + h12 y + h13 − u h31 x − u h32 y = u
   *   h21 x + h22 y + h23 − v h31 x − v h32 y = v
   */
  for (i = 0; i < 4; i++)
    {
      double x = source->p[i].x;
      double y = source->p[i].y;
      double u = dest->p[i].x;
      double v = dest->p[i].y;
      double *r0 = a[2 * i];
      double *r1 = a[2 * i + 1];

      r0[0] = x; r0[1] = y; r0[2] = 1; r0[3] = 0; r0[4] = 0; r0[5] = 0;
      r0[6] = -u * x; r0[7] = -u * y; r0[8] = u;

      r1[0] = 0; r1[1] = 0; r1[2] = 0; r1[3] = x; r1[4] = y; r1[5] = 1;
      r1[6] = -v * x; r1[7] = -v * y; r1[8] = v;
    }

  /* Forward elimination */
  for (col = 0; col < 8; col++)
    {
      int pivot = col;

      for (row = col + 1; row < 8; row++)
        if (fabs (a[row][col]) > fabs (a[pivot][col]))
          pivot = row;

      if (fabs (a[pivot][col]) < 1e-12)
        return FALSE;

      if (pivot != col)
        {
          double tmp[9];

          memcpy (tmp, a[col], sizeof tmp);
          memcpy (a[col], a[pivot], sizeof tmp);
          memcpy (a[pivot], tmp, sizeof tmp);
        }

      for (row = col + 1; row < 8; row++)
        {
          double factor = a[row][col] / a[col][col];
          int k;

          for (k = col; k < 9; k++)
            a[row][k] -= factor * a[col][k];
        }
    }

  /* Back substitution */
  for (row = 7; row >= 0; row--)
    {
      double sum = a[row][8];

      for (col = row + 1; col < 8; col++)
        sum -= a[row][col] * h[col];

      h[row] = sum / a[row][row];
    }

  matrix->coeff[0][0] = h[0];
  matrix->coeff[0][1] = h[1];
  matrix->coeff[0][2] = h[2];
  matrix->coeff[1][0] = h[3];
  matrix->coeff[1][1] = h[4];
  matrix->coeff[1][2] = h[5];
  matrix->coeff[2][0] = h[6];
  matrix->coeff[2][1] = h[7];
  matrix->coeff[2][2] = 1.0;

  return TRUE;
}

static int
bm_point_compare (const void *a,
                  const void *b)
{
  const BmPoint *p = a;
  const BmPoint *q = b;

  if (p->x != q->x)
    return p->x < q->x ? -1 : 1;
  if (p->y != q->y)
    return p->y < q->y ? -1 : 1;

  return 0;
}

static double
bm_point_cross (const BmPoint *o,
                const BmPoint *a,
                const BmPoint *b)
{
  return (a->x - o->x) * (b->y - o->y) - (a->y - o->y) * (b->x - o->x);
}

/**
 * bm_polygon_convex_hull:
 * @points: (array length=n_points): points to wrap; sorted in place
 * @n_points: number of points, at most %BM_MAX_POLYGON_POINTS
 * @hull: (out caller-allocates): return location for the hull
 *
 * Computes the convex hull with Andrew's monotone chain.
 */
static void
bm_polygon_convex_hull (BmPoint   *points,
                        int        n_points,
                        BmPolygon *hull)
{
  BmPoint chain[2 * BM_MAX_POLYGON_POINTS];
  int n = 0;
  int lower_end;
  int i;

  g_return_if_fail (n_points <= BM_MAX_POLYGON_POINTS);

  qsort (points, n_points, sizeof (BmPoint), bm_point_compare);

  for (i = 0; i < n_points; i++)
    {
      while (n >= 2 && bm_point_cross (&chain[n - 2], &chain[n - 1],
                                       &points[i]) <= 0)
        n--;
      chain[n++] = points[i];
    }

  lower_end = n + 1;
  for (i = n_points - 2; i >= 0; i--)
    {
      while (n >= lower_end && bm_point_cross (&chain[n - 2], &chain[n - 1],
                                               &points[i]) <= 0)
        n--;
      chain[n++] = points[i];
    }

  /* The last point repeats the first */
  hull->n_points = MAX (n - 1, 0);
  memcpy (hull->points, chain, hull->n_points * sizeof (BmPoint));
}


/* Stratum 2: Scene
 *
 * A pinhole camera looking down +z, flat textured faces in 3D, an upright
 * box built from them, and the box's shadow on the floor. Depends on
 * stratum 1 only.
 */

#define BM_MAX_FACES    6
#define BM_MAX_BANDS    32
#define BM_MAX_TEXTURES 4

typedef struct _BmCamera  BmCamera;
typedef struct _BmShade   BmShade;
typedef struct _BmTexture BmTexture;
typedef struct _BmFace    BmFace;
typedef struct _BmLight  BmLight;
typedef struct _BmScene  BmScene;
typedef struct _BmBand   BmBand;

struct _BmCamera
{
  double focal;             /* pixels; larger is flatter                */
  double offset_x;          /* screen position of x = 0                 */
  double eye_y;             /* screen position of y = 0                 */
};

/* Linear brightness ramp across a face, from its left edge to its right
 * edge in source space. Equal values give a flat tint.
 */
struct _BmShade
{
  double start;             /* 0..1 multiplier at the left edge         */
  double end;               /* 0..1 multiplier at the right edge        */
};

typedef enum
{
  BM_TEXTURE_NOISE,         /* generated gray noise, optionally streaked */
  BM_TEXTURE_IMAGE,         /* a texture file, stretched to the face     */
  BM_TEXTURE_CREASE,        /* a vertical groove, such as a hinge score  */
} BmTextureKind;

/* A detail map multiplied onto a face. The map is read as gray and
 * centered on mid-gray, so the face is multiplied by 1 ± @strength:
 * texture adds detail without changing overall brightness. Layers apply
 * in order, before shading.
 */
struct _BmTexture
{
  BmTextureKind kind;
  const char   *image;      /* IMAGE: texture file                       */
  double        strength;   /* 0..1 swing of the multiplier              */
  double        streak;     /* NOISE: blur length in source units;
                             * 0 gives isotropic grain                   */
  double        angle;      /* NOISE: streak direction in source space,
                             * degrees; 0 runs left to right             */
  guint         seed;       /* NOISE: fixed for reproducible output      */
  double        offset;     /* CREASE: center, from the left edge, in
                             * source units                              */
  double        width;      /* CREASE: groove width, in source units     */
  gboolean      flip;       /* CREASE: lit wall on the right; set by
                             * bm_scene_light()                          */
};

/* A flat, four-cornered face. When @image is %NULL the face is filled with
 * @color instead, for surfaces without artwork such as a page block.
 */
struct _BmFace
{
  const char *image;        /* artwork file, or NULL                    */
  double      color[3];     /* linear RGB, used when image is NULL      */
  double      width;        /* source size, in scene units              */
  double      height;
  BmVec3      corners[4];   /* TL, TR, BR, BL in scene space            */
  BmShade     shade;

  BmTexture   textures[BM_MAX_TEXTURES];
  int         n_textures;
};

/* A directional light from the front, above the scene. Its azimuth and
 * shadow length fix both the floor shadow and the light's elevation
 * (a shadow L times the object height long means cot(elevation) = L), so
 * faces are lit consistently with the shadow they cast.
 */
struct _BmLight
{
  double azimuth;           /* degrees from straight back, + is left    */
  double length;            /* shadow length, fraction of object height */
  double ambient;           /* brightness of faces turned from the light */
  double diffuse;           /* added brightness for faces facing it     */
};

/* Sides of a box, named as seen with its front toward the camera */
typedef enum
{
  BM_BOX_FRONT,
  BM_BOX_BACK,
  BM_BOX_LEFT,
  BM_BOX_RIGHT,
  BM_BOX_TOP,
  BM_BOX_BOTTOM,
  BM_BOX_N_SIDES,
} BmBoxSide;

struct _BmScene
{
  BmCamera  camera;
  BmFace    faces[BM_MAX_FACES];
  int       n_faces;

  BmPolygon footprint;      /* floor contact; points are (x, z)         */
  double    object_height;
  double    floor_y;        /* y of the floor; > 0 is below the eye     */
};

struct _BmBand
{
  BmPolygon polygon;        /* screen space                             */
  double    level;          /* 0..1 mask intensity                      */
};

/**
 * bm_camera_project:
 * @camera: a `BmCamera`
 * @point: a scene-space point
 *
 * Projects @point with k = f / (f + z).
 *
 * Returns: the screen-space point
 */
static BmPoint
bm_camera_project (const BmCamera *camera,
                   const BmVec3   *point)
{
  double k = camera->focal / (camera->focal + point->z);

  return (BmPoint) { point->x * k + camera->offset_x,
                     point->y * k + camera->eye_y };
}

/**
 * bm_face_project:
 * @face: a `BmFace`
 * @camera: a `BmCamera`
 *
 * Projects the corners of @face.
 *
 * Returns: the face's screen-space quad
 */
static BmQuad
bm_face_project (const BmFace   *face,
                 const BmCamera *camera)
{
  BmQuad quad;
  int i;

  for (i = 0; i < 4; i++)
    quad.p[i] = bm_camera_project (camera, &face->corners[i]);

  return quad;
}

/**
 * bm_face_add_texture:
 * @face: a `BmFace`
 * @texture: the texture layer to append
 *
 * Appends a texture layer to @face; layers apply in the order added.
 *
 * Returns: %FALSE if @face already has %BM_MAX_TEXTURES layers
 */
static gboolean
bm_face_add_texture (BmFace          *face,
                     const BmTexture *texture)
{
  if (face->n_textures >= BM_MAX_TEXTURES)
    return FALSE;

  face->textures[face->n_textures++] = *texture;

  return TRUE;
}

/**
 * bm_face_is_visible:
 * @face: a `BmFace`
 * @camera: a `BmCamera`
 *
 * Tests whether @face points toward the camera. Corners are listed
 * clockwise as seen from the front, so a visible face projects with a
 * positive signed area (y grows down) and a face seen from behind
 * projects with a negative one.
 *
 * Returns: %TRUE if the face is visible
 */
static gboolean
bm_face_is_visible (const BmFace   *face,
                    const BmCamera *camera)
{
  BmQuad quad = bm_face_project (face, camera);
  double area = 0.0;
  int i;

  for (i = 0; i < 4; i++)
    {
      const BmPoint *a = &quad.p[i];
      const BmPoint *b = &quad.p[(i + 1) % 4];

      area += a->x * b->y - b->x * a->y;
    }

  return area > 0.0;
}

/**
 * bm_scene_add_box:
 * @scene: a `BmScene`
 * @sides: (array fixed-size=6): a face per `BmBoxSide`, each nullable;
 *   only image, color, shade and textures are read
 * @width: box width, along the front (x)
 * @depth: box depth, front to back (z)
 * @height: box height (y)
 * @angle: rotation about the box's vertical center line, in degrees
 * @eye_height: height of the camera above the floor, in scene units
 *
 * Adds a box standing on the floor. At @angle 0 the front faces the
 * camera; positive angles swing the left side toward the camera, 90 shows
 * the left side square on, and 180 shows the back. The box is moved back
 * so that its nearest point sits at z = 0, where the projection scale is
 * exactly 1.
 *
 * Each face's source image is upright as seen from outside the box. The
 * top shows the front edge at the bottom of its image; the bottom shows
 * it at the top. Fills in the face corners and sizes, the footprint, the
 * object height and the floor position. Faces turned away from the
 * camera are still added; check them with bm_face_is_visible().
 */
static void
bm_scene_add_box (BmScene       *scene,
                  const BmFace  *sides[BM_BOX_N_SIDES],
                  double         width,
                  double         depth,
                  double         height,
                  double         angle,
                  double         eye_height)
{
  const double w = width;
  const double d = depth;
  const double y0 = eye_height - height;    /* top    */
  const double y1 = eye_height;             /* bottom */
  /* Local corners per side: TL, TR, BR, BL as (x, y, z), x from the left
   * side, z from the front.
   */
  const BmVec3 local[BM_BOX_N_SIDES][4] = {
    [BM_BOX_FRONT]  = { { 0, y0, 0 }, { w, y0, 0 },
                        { w, y1, 0 }, { 0, y1, 0 } },
    [BM_BOX_BACK]   = { { w, y0, d }, { 0, y0, d },
                        { 0, y1, d }, { w, y1, d } },
    [BM_BOX_LEFT]   = { { 0, y0, d }, { 0, y0, 0 },
                        { 0, y1, 0 }, { 0, y1, d } },
    [BM_BOX_RIGHT]  = { { w, y0, 0 }, { w, y0, d },
                        { w, y1, d }, { w, y1, 0 } },
    [BM_BOX_TOP]    = { { 0, y0, d }, { w, y0, d },
                        { w, y0, 0 }, { 0, y0, 0 } },
    [BM_BOX_BOTTOM] = { { 0, y1, 0 }, { w, y1, 0 },
                        { w, y1, d }, { 0, y1, d } },
  };
  const double size[BM_BOX_N_SIDES][2] = {
    [BM_BOX_FRONT]  = { w, height },
    [BM_BOX_BACK]   = { w, height },
    [BM_BOX_LEFT]   = { d, height },
    [BM_BOX_RIGHT]  = { d, height },
    [BM_BOX_TOP]    = { w, d },
    [BM_BOX_BOTTOM] = { w, d },
  };
  const BmPoint floor[4] = { { 0, 0 }, { w, 0 }, { w, d }, { 0, d } };
  double t = BM_DEG_TO_RAD (angle);
  double c = cos (t);
  double s = sin (t);
  double z_shift = G_MAXDOUBLE;
  int side, i;

  /* Rotation about the vertical line through the box's center:
   *   x' = (x − w/2) cos t − (z − d/2) sin t
   *   z' = (x − w/2) sin t + (z − d/2) cos t
   * then a z shift that puts the nearest corner at z = 0.
   */
#define BM_ROTATE_X(x, z) (((x) - w / 2.0) * c - ((z) - d / 2.0) * s)
#define BM_ROTATE_Z(x, z) (((x) - w / 2.0) * s + ((z) - d / 2.0) * c)

  for (i = 0; i < 4; i++)
    z_shift = MIN (z_shift, BM_ROTATE_Z (floor[i].x, floor[i].y));
  z_shift = -z_shift;

  for (side = 0; side < BM_BOX_N_SIDES; side++)
    {
      BmFace *f;

      if (sides[side] == NULL)
        continue;

      g_return_if_fail (scene->n_faces < BM_MAX_FACES);

      f = &scene->faces[scene->n_faces++];
      *f = *sides[side];
      f->width = size[side][0];
      f->height = size[side][1];

      for (i = 0; i < 4; i++)
        {
          const BmVec3 *p = &local[side][i];

          f->corners[i] = (BmVec3) { BM_ROTATE_X (p->x, p->z),
                                     p->y,
                                     BM_ROTATE_Z (p->x, p->z) + z_shift };
        }
    }

  scene->footprint.n_points = 4;
  for (i = 0; i < 4; i++)
    scene->footprint.points[i] =
      (BmPoint) { BM_ROTATE_X (floor[i].x, floor[i].y),
                  BM_ROTATE_Z (floor[i].x, floor[i].y) + z_shift };

#undef BM_ROTATE_X
#undef BM_ROTATE_Z

  scene->object_height = height;
  scene->floor_y = y1;
}

static BmVec3
bm_vec3_sub (const BmVec3 *a,
             const BmVec3 *b)
{
  return (BmVec3) { a->x - b->x, a->y - b->y, a->z - b->z };
}

static BmVec3
bm_vec3_cross (const BmVec3 *a,
               const BmVec3 *b)
{
  return (BmVec3) { a->y * b->z - a->z * b->y,
                    a->z * b->x - a->x * b->z,
                    a->x * b->y - a->y * b->x };
}

static double
bm_vec3_dot (const BmVec3 *a,
             const BmVec3 *b)
{
  return a->x * b->x + a->y * b->y + a->z * b->z;
}

static BmVec3
bm_vec3_normalize (const BmVec3 *v)
{
  double length = sqrt (bm_vec3_dot (v, v));

  if (length == 0.0)
    return *v;

  return (BmVec3) { v->x / length, v->y / length, v->z / length };
}

/**
 * bm_light_direction:
 * @light: a `BmLight`
 *
 * Returns: the unit vector pointing toward the light
 */
static BmVec3
bm_light_direction (const BmLight *light)
{
  double azimuth = BM_DEG_TO_RAD (light->azimuth);
  double elevation = atan2 (1.0, MAX (light->length, 1e-6));
  BmVec3 v = { sin (azimuth) * cos (elevation),
               -sin (elevation),
               -cos (azimuth) * cos (elevation) };

  return bm_vec3_normalize (&v);
}

/**
 * bm_scene_light:
 * @scene: a `BmScene`
 * @light: the light
 * @falloff: 0..1 brightness at the scene's far depth, relative to its
 *   nearest point
 *
 * Lights every face: a Lambert term from the face's outward normal (its
 * corners are clockwise seen from outside, so the normal is
 * (BL − TL) × (TR − TL)), times a depth falloff measured at the face's
 * left and right edges. The result multiplies each face's own shade ramp.
 * Creases are flipped so their lit wall faces the light.
 */
static void
bm_scene_light (BmScene       *scene,
                const BmLight *light,
                double         falloff)
{
  BmVec3 toward = bm_light_direction (light);
  double z_min = G_MAXDOUBLE;
  double z_max = -G_MAXDOUBLE;
  int i, j;

  for (i = 0; i < scene->n_faces; i++)
    for (j = 0; j < 4; j++)
      {
        z_min = MIN (z_min, scene->faces[i].corners[j].z);
        z_max = MAX (z_max, scene->faces[i].corners[j].z);
      }

  for (i = 0; i < scene->n_faces; i++)
    {
      BmFace *f = &scene->faces[i];
      BmVec3 right = bm_vec3_sub (&f->corners[1], &f->corners[0]);
      BmVec3 down = bm_vec3_sub (&f->corners[3], &f->corners[0]);
      BmVec3 cross = bm_vec3_cross (&down, &right);
      BmVec3 normal = bm_vec3_normalize (&cross);
      double lit = light->ambient +
                   light->diffuse * MAX (bm_vec3_dot (&normal, &toward), 0.0);
      double z_left = (f->corners[0].z + f->corners[3].z) / 2.0;
      double z_right = (f->corners[1].z + f->corners[2].z) / 2.0;
      double span = MAX (z_max - z_min, 1e-6);

      lit = MIN (lit, 1.0);
      f->shade.start *= lit * (1.0 - (1.0 - falloff) * (z_left - z_min) / span);
      f->shade.end *= lit * (1.0 - (1.0 - falloff) * (z_right - z_min) / span);

      for (j = 0; j < f->n_textures; j++)
        if (f->textures[j].kind == BM_TEXTURE_CREASE)
          f->textures[j].flip = bm_vec3_dot (&right, &toward) < 0.0;
    }
}

/**
 * bm_scene_center:
 * @scene: a `BmScene`
 * @canvas_width: canvas width
 * @canvas_height: canvas height
 *
 * Sets the camera offsets so the visible faces are centered on the
 * canvas in both directions.
 */
static void
bm_scene_center (BmScene *scene,
                 int      canvas_width,
                 int      canvas_height)
{
  double min_x = G_MAXDOUBLE;
  double max_x = -G_MAXDOUBLE;
  double min_y = G_MAXDOUBLE;
  double max_y = -G_MAXDOUBLE;
  gboolean found = FALSE;
  int i, j;

  scene->camera.offset_x = 0.0;
  scene->camera.eye_y = 0.0;

  for (i = 0; i < scene->n_faces; i++)
    {
      BmQuad quad;

      if (!bm_face_is_visible (&scene->faces[i], &scene->camera))
        continue;

      quad = bm_face_project (&scene->faces[i], &scene->camera);
      found = TRUE;

      for (j = 0; j < 4; j++)
        {
          min_x = MIN (min_x, quad.p[j].x);
          max_x = MAX (max_x, quad.p[j].x);
          min_y = MIN (min_y, quad.p[j].y);
          max_y = MAX (max_y, quad.p[j].y);
        }
    }

  if (!found)
    return;

  scene->camera.offset_x = (canvas_width - (max_x - min_x)) / 2.0 - min_x;
  scene->camera.eye_y = (canvas_height - (max_y - min_y)) / 2.0 - min_y;
}

/**
 * bm_scene_shadow_bands:
 * @scene: a `BmScene`
 * @light: the light
 * @t_max: fraction of the full shadow length for the outermost band
 * @n_steps: number of bands, at most %BM_MAX_BANDS
 * @opacity: mask level of the innermost band
 * @bands: (out caller-allocates) (array): return location for the bands
 *
 * Builds floor-shadow polygons: the hull of the footprint and the
 * footprint offset along the light, for t = @t_max … @t_max / N, projected
 * to the screen. Bands come largest first, each brighter than the last,
 * so that painting them in order leaves the base densest.
 *
 * Returns: the number of bands written
 */
static int
bm_scene_shadow_bands (const BmScene *scene,
                       const BmLight *light,
                       double         t_max,
                       int            n_steps,
                       double         opacity,
                       BmBand        *bands)
{
  const BmPolygon *foot = &scene->footprint;
  double azimuth = BM_DEG_TO_RAD (light->azimuth);
  double length = light->length * scene->object_height;
  double dx = -sin (azimuth) * length;
  double dz = cos (azimuth) * length;
  int n_bands = 0;
  int step, i;

  g_return_val_if_fail (n_steps <= BM_MAX_BANDS, 0);
  g_return_val_if_fail (2 * foot->n_points <= BM_MAX_POLYGON_POINTS, 0);

  for (step = n_steps; step >= 1; step--)
    {
      double t = t_max * step / n_steps;
      BmPoint points[BM_MAX_POLYGON_POINTS];
      BmPolygon hull;
      BmBand *band = &bands[n_bands++];

      for (i = 0; i < foot->n_points; i++)
        {
          points[i] = foot->points[i];
          points[i + foot->n_points] =
            (BmPoint) { foot->points[i].x + t * dx,
                        foot->points[i].y + t * dz };
        }

      bm_polygon_convex_hull (points, 2 * foot->n_points, &hull);

      band->polygon.n_points = hull.n_points;
      for (i = 0; i < hull.n_points; i++)
        {
          BmVec3 floor_point = { hull.points[i].x, scene->floor_y,
                                 hull.points[i].y };

          band->polygon.points[i] = bm_camera_project (&scene->camera,
                                                       &floor_point);
        }

      band->level = opacity * pow ((double) (n_steps - step + 1) / n_steps,
                                   1.6);
    }

  return n_bands;
}


/* Stratum 3: Render
 *
 * Generic GEGL graph builders. Every helper adds nodes under @root and
 * returns the node whose "output" pad carries the result; `gegl:over`
 * composites its "aux" pad over its "input" pad. Depends on strata 1
 * and 2.
 */

/**
 * bm_render_gray:
 * @level: 0..1 linear intensity
 * @alpha: 0..1 alpha
 *
 * Returns: (transfer full): a new `GeglColor`
 */
static GeglColor *
bm_render_gray (double level,
                double alpha)
{
  GeglColor *color = gegl_color_new (NULL);

  gegl_color_set_rgba (color, level, level, level, alpha);

  return color;
}

/**
 * bm_render_load:
 * @root: the graph root
 * @path: image file to load
 * @width: target width
 * @height: target height
 *
 * Adds `gegl:load` followed by `gegl:scale-size`.
 *
 * Returns: (transfer none): the scale node
 */
static GeglNode *
bm_render_load (GeglNode   *root,
                const char *path,
                double      width,
                double      height)
{
  GeglNode *load;
  GeglNode *scale;

  load = gegl_node_new_child (root,
                              "operation", "gegl:load",
                              "path", path,
                              NULL);
  scale = gegl_node_new_child (root,
                               "operation", "gegl:scale-size",
                               "x", width,
                               "y", height,
                               NULL);
  gegl_node_link (load, scale);

  return scale;
}

/**
 * bm_render_solid:
 * @root: the graph root
 * @rgb: (array fixed-size=3): linear RGB color
 * @width: rectangle width
 * @height: rectangle height
 *
 * Adds a solid-color rectangle at the origin (`gegl:color` cropped).
 *
 * Returns: (transfer none): the crop node
 */
static GeglNode *
bm_render_solid (GeglNode     *root,
                 const double *rgb,
                 double        width,
                 double        height)
{
  GeglColor *color = gegl_color_new (NULL);
  GeglNode *fill;
  GeglNode *crop;

  gegl_color_set_rgba (color, rgb[0], rgb[1], rgb[2], 1.0);
  fill = gegl_node_new_child (root,
                              "operation", "gegl:color",
                              "value", color,
                              NULL);
  g_object_unref (color);

  crop = gegl_node_new_child (root,
                              "operation", "gegl:crop",
                              "x", 0.0,
                              "y", 0.0,
                              "width", width,
                              "height", height,
                              NULL);
  gegl_node_link (fill, crop);

  return crop;
}

/**
 * bm_render_shade:
 * @root: the graph root
 * @input: the node to shade
 * @shade: the brightness ramp
 * @width: source width
 * @height: source height
 *
 * Multiplies @input by a horizontal gray ramp (`gegl:linear-gradient`,
 * cropped to the face, as the aux of `gegl:multiply`).
 *
 * Returns: (transfer none): the multiply node
 */
static GeglNode *
bm_render_shade (GeglNode      *root,
                 GeglNode      *input,
                 const BmShade *shade,
                 double         width,
                 double         height)
{
  GeglColor *start = bm_render_gray (shade->start, 1.0);
  GeglColor *end = bm_render_gray (shade->end, 1.0);
  GeglNode *gradient;
  GeglNode *crop;
  GeglNode *multiply;

  gradient = gegl_node_new_child (root,
                                  "operation", "gegl:linear-gradient",
                                  "start-x", 0.0,
                                  "start-y", 0.0,
                                  "end-x", width,
                                  "end-y", 0.0,
                                  "start-color", start,
                                  "end-color", end,
                                  NULL);
  g_object_unref (start);
  g_object_unref (end);

  crop = gegl_node_new_child (root,
                              "operation", "gegl:crop",
                              "x", 0.0,
                              "y", 0.0,
                              "width", width,
                              "height", height,
                              NULL);
  multiply = gegl_node_new_child (root,
                                  "operation", "gegl:multiply",
                                  NULL);

  gegl_node_link (gradient, crop);
  gegl_node_link (input, multiply);
  gegl_node_connect (multiply, "aux", crop, "output");

  return multiply;
}

/**
 * bm_render_noise:
 * @root: the graph root
 * @texture: a %BM_TEXTURE_NOISE texture
 * @width: map width
 * @height: map height
 *
 * Generates a gray noise map stretched to the full 0..1 range. With a
 * nonzero streak the noise is motion-blurred along the texture's angle,
 * which turns grain into fine parallel lines, such as page edges. The
 * noise is generated with a margin so the blur does not fade at the
 * edges.
 *
 * Returns: (transfer none): the map node
 */
static GeglNode *
bm_render_noise (GeglNode        *root,
                 const BmTexture *texture,
                 double           width,
                 double           height)
{
  GeglColor *gray = bm_render_gray (0.5, 1.0);
  double margin = ceil (texture->streak) + 2.0;
  GeglNode *last;
  GeglNode *node;

  node = gegl_node_new_child (root,
                              "operation", "gegl:color",
                              "value", gray,
                              NULL);
  g_object_unref (gray);

  last = gegl_node_new_child (root,
                              "operation", "gegl:crop",
                              "x", -margin,
                              "y", -margin,
                              "width", width + 2.0 * margin,
                              "height", height + 2.0 * margin,
                              NULL);
  gegl_node_link (node, last);

  node = gegl_node_new_child (root,
                              "operation", "gegl:noise-rgb",
                              "independent", FALSE,
                              "red", 0.5,
                              "seed", texture->seed,
                              NULL);
  gegl_node_link (last, node);
  last = node;

  if (texture->streak > 0.0)
    {
      node = gegl_node_new_child (root,
                                  "operation", "gegl:motion-blur-linear",
                                  "length", texture->streak,
                                  "angle", texture->angle,
                                  NULL);
      gegl_node_link (last, node);
      last = node;
    }

  node = gegl_node_new_child (root,
                              "operation", "gegl:crop",
                              "x", 0.0,
                              "y", 0.0,
                              "width", width,
                              "height", height,
                              NULL);
  gegl_node_link (last, node);
  last = node;

  node = gegl_node_new_child (root,
                              "operation", "gegl:stretch-contrast",
                              NULL);
  gegl_node_link (last, node);

  return node;
}

/**
 * bm_render_rect_path:
 * @x: left edge
 * @y: top edge
 * @width: rectangle width
 * @height: rectangle height
 *
 * Returns: (transfer full): a new closed rectangular `GeglPath`
 */
static GeglPath *
bm_render_rect_path (double x,
                     double y,
                     double width,
                     double height)
{
  char x0[G_ASCII_DTOSTR_BUF_SIZE];
  char y0[G_ASCII_DTOSTR_BUF_SIZE];
  char x1[G_ASCII_DTOSTR_BUF_SIZE];
  char y1[G_ASCII_DTOSTR_BUF_SIZE];
  g_autofree char *d = NULL;

  g_ascii_dtostr (x0, sizeof x0, x);
  g_ascii_dtostr (y0, sizeof y0, y);
  g_ascii_dtostr (x1, sizeof x1, x + width);
  g_ascii_dtostr (y1, sizeof y1, y + height);
  d = g_strdup_printf ("M %s %s L %s %s L %s %s L %s %s Z",
                       x0, y0, x1, y0, x1, y1, x0, y1);

  return gegl_path_new_from_string (d);
}

/**
 * bm_render_crease:
 * @root: the graph root
 * @texture: a %BM_TEXTURE_CREASE texture
 * @width: map width
 * @height: map height
 *
 * Generates a groove profile on mid-gray: a lit wall (white) left of the
 * crease line and a shaded wall (black) right of it, for light from the
 * face's right, or the reverse when the texture is flipped. A horizontal
 * blur rounds the profile.
 *
 * Returns: (transfer none): the map node
 */
static GeglNode *
bm_render_crease (GeglNode        *root,
                  const BmTexture *texture,
                  double           width,
                  double           height)
{
  double half = texture->width / 2.0;
  const struct
  {
    double x;
    double level;
  } walls[] = {
    { texture->offset - half, texture->flip ? 0.0 : 1.0 },
    { texture->offset,        texture->flip ? 1.0 : 0.0 },
  };
  GeglNode *last;
  GeglNode *node;
  gsize i;

  last = bm_render_solid (root, (double[]) { 0.5, 0.5, 0.5 }, width, height);

  for (i = 0; i < G_N_ELEMENTS (walls); i++)
    {
      GeglColor *color = bm_render_gray (walls[i].level, 1.0);
      GeglPath *path = bm_render_rect_path (walls[i].x, 0.0, half, height);

      node = gegl_node_new_child (root,
                                  "operation", "gegl:fill-path",
                                  "color", color,
                                  "d", path,
                                  NULL);
      g_object_unref (color);
      g_object_unref (path);

      gegl_node_link (last, node);
      last = node;
    }

  node = gegl_node_new_child (root,
                              "operation", "gegl:gaussian-blur",
                              "std-dev-x", MAX (half / 2.0, 0.5),
                              "std-dev-y", 0.0,
                              NULL);
  gegl_node_link (last, node);

  return node;
}

/**
 * bm_render_texture:
 * @root: the graph root
 * @input: the node to texture
 * @texture: the texture layer
 * @width: face width
 * @height: face height
 *
 * Builds the texture's gray map, centers it on mid-gray as a multiplier
 * of 1 ± strength, and multiplies @input by it. The math operations
 * leave @input's alpha untouched.
 *
 * Returns: (transfer none): the multiply node
 */
static GeglNode *
bm_render_texture (GeglNode        *root,
                   GeglNode        *input,
                   const BmTexture *texture,
                   double           width,
                   double           height)
{
  GeglNode *map;
  GeglNode *center;
  GeglNode *scale;
  GeglNode *offset;
  GeglNode *multiply;

  switch (texture->kind)
    {
    case BM_TEXTURE_IMAGE:
      {
        GeglNode *gray = gegl_node_new_child (root,
                                              "operation", "gegl:gray",
                                              NULL);

        map = bm_render_load (root, texture->image, width, height);
        gegl_node_link (map, gray);
        map = gray;
      }
      break;

    case BM_TEXTURE_CREASE:
      map = bm_render_crease (root, texture, width, height);
      break;

    case BM_TEXTURE_NOISE:
    default:
      map = bm_render_noise (root, texture, width, height);
      break;
    }

  /* m → 1 + 2 · strength · (m − 0.5) */
  center = gegl_node_new_child (root,
                                "operation", "gegl:subtract",
                                "value", 0.5,
                                NULL);
  scale = gegl_node_new_child (root,
                               "operation", "gegl:multiply",
                               "value", 2.0 * texture->strength,
                               NULL);
  offset = gegl_node_new_child (root,
                                "operation", "gegl:add",
                                "value", 1.0,
                                NULL);
  multiply = gegl_node_new_child (root,
                                  "operation", "gegl:multiply",
                                  NULL);

  gegl_node_link_many (map, center, scale, offset, NULL);
  gegl_node_link (input, multiply);
  gegl_node_connect (multiply, "aux", offset, "output");

  return multiply;
}

/**
 * bm_render_warp:
 * @root: the graph root
 * @input: the node to transform
 * @matrix: a projective matrix
 *
 * Adds a `gegl:transform` node. Its "transform" property accepts all nine
 * coefficients as `matrix(…)`; always produce that string with
 * gegl_matrix3_to_string(), since the string lists them column by column.
 *
 * Returns: (transfer none): the transform node
 */
static GeglNode *
bm_render_warp (GeglNode          *root,
                GeglNode          *input,
                const GeglMatrix3 *matrix)
{
  g_autofree char *transform = gegl_matrix3_to_string (matrix);
  GeglNode *node;

  node = gegl_node_new_child (root,
                              "operation", "gegl:transform",
                              "transform", transform,
                              NULL);
  gegl_node_link (input, node);

  return node;
}

/**
 * bm_render_band_path:
 * @band: a `BmBand`
 *
 * Converts a band into an SVG-style "M … L … Z" path.
 *
 * Returns: (transfer full): a new `GeglPath`
 */
static GeglPath *
bm_render_band_path (const BmBand *band)
{
  g_autoptr (GString) d = g_string_new (NULL);
  char x[G_ASCII_DTOSTR_BUF_SIZE];
  char y[G_ASCII_DTOSTR_BUF_SIZE];
  int i;

  for (i = 0; i < band->polygon.n_points; i++)
    {
      g_ascii_dtostr (x, sizeof x, band->polygon.points[i].x);
      g_ascii_dtostr (y, sizeof y, band->polygon.points[i].y);
      g_string_append_printf (d, "%c %s %s ", i == 0 ? 'M' : 'L', x, y);
    }
  g_string_append (d, "Z");

  return gegl_path_new_from_string (d->str);
}

/**
 * bm_render_mask:
 * @root: the graph root
 * @bands: (array length=n_bands): bands to paint, in order
 * @n_bands: number of bands
 * @blur: Gaussian standard deviation, in pixels
 * @width: canvas width
 * @height: canvas height
 *
 * Paints polygons onto black with chained `gegl:fill-path` nodes, each in
 * gray at its band's level, then blurs the result.
 *
 * Returns: (transfer none): the blurred mask node
 */
static GeglNode *
bm_render_mask (GeglNode     *root,
                const BmBand *bands,
                int           n_bands,
                double        blur,
                int           width,
                int           height)
{
  GeglColor *black = bm_render_gray (0.0, 1.0);
  GeglNode *last;
  GeglNode *node;
  int i;

  node = gegl_node_new_child (root,
                              "operation", "gegl:color",
                              "value", black,
                              NULL);
  g_object_unref (black);
  last = gegl_node_new_child (root,
                              "operation", "gegl:crop",
                              "x", 0.0,
                              "y", 0.0,
                              "width", (double) width,
                              "height", (double) height,
                              NULL);
  gegl_node_link (node, last);

  for (i = 0; i < n_bands; i++)
    {
      GeglColor *gray = bm_render_gray (bands[i].level, 1.0);
      GeglPath *path = bm_render_band_path (&bands[i]);

      node = gegl_node_new_child (root,
                                  "operation", "gegl:fill-path",
                                  "color", gray,
                                  "d", path,
                                  NULL);
      g_object_unref (gray);
      g_object_unref (path);
      gegl_node_link (last, node);
      last = node;
    }

  node = gegl_node_new_child (root,
                              "operation", "gegl:gaussian-blur",
                              "std-dev-x", blur,
                              "std-dev-y", blur,
                              NULL);
  gegl_node_link (last, node);

  return node;
}

/**
 * bm_render_shadow:
 * @root: the graph root
 * @masks: (array length=n_masks): mask nodes to merge
 * @n_masks: number of masks, at least one
 * @width: canvas width
 * @height: canvas height
 *
 * Merges masks with `svg:lighten` and returns a black layer whose opacity
 * follows the result (`gegl:opacity` with the mask as aux).
 *
 * The mask only ever reaches the alpha channel: the layer's color stays
 * exactly #000000. Rendered over a transparent background, the shadow is
 * therefore pure black with varying alpha, never gray, and composites
 * cleanly onto dark backgrounds.
 *
 * Returns: (transfer none): the shadow layer node
 */
static GeglNode *
bm_render_shadow (GeglNode  *root,
                  GeglNode **masks,
                  int        n_masks,
                  int        width,
                  int        height)
{
  GeglColor *black = bm_render_gray (0.0, 1.0);
  GeglNode *mask;
  GeglNode *color;
  GeglNode *crop;
  GeglNode *opacity;
  int i;

  g_return_val_if_fail (n_masks >= 1, NULL);

  mask = masks[0];
  for (i = 1; i < n_masks; i++)
    {
      GeglNode *lighten = gegl_node_new_child (root,
                                               "operation", "svg:lighten",
                                               NULL);

      gegl_node_link (mask, lighten);
      gegl_node_connect (lighten, "aux", masks[i], "output");
      mask = lighten;
    }

  color = gegl_node_new_child (root,
                               "operation", "gegl:color",
                               "value", black,
                               NULL);
  g_object_unref (black);
  crop = gegl_node_new_child (root,
                              "operation", "gegl:crop",
                              "x", 0.0,
                              "y", 0.0,
                              "width", (double) width,
                              "height", (double) height,
                              NULL);
  opacity = gegl_node_new_child (root,
                                 "operation", "gegl:opacity",
                                 NULL);

  gegl_node_link_many (color, crop, opacity, NULL);
  gegl_node_connect (opacity, "aux", mask, "output");

  return opacity;
}

/**
 * bm_render_face:
 * @root: the graph root
 * @face: the face to render
 * @camera: the camera
 *
 * Loads (or fills), textures, shades and warps @face into place on the
 * canvas.
 *
 * Returns: (transfer none) (nullable): the face's output node, or %NULL if
 *   the face is degenerate
 */
static GeglNode *
bm_render_face (GeglNode       *root,
                const BmFace   *face,
                const BmCamera *camera)
{
  BmQuad source = bm_quad_from_rect (face->width, face->height);
  BmQuad dest = bm_face_project (face, camera);
  GeglMatrix3 matrix;
  GeglNode *node;
  int i;

  if (!bm_homography_from_quads (&source, &dest, &matrix))
    return NULL;

  if (face->image != NULL)
    node = bm_render_load (root, face->image, face->width, face->height);
  else
    node = bm_render_solid (root, face->color, face->width, face->height);

  for (i = 0; i < face->n_textures; i++)
    node = bm_render_texture (root, node, &face->textures[i],
                              face->width, face->height);

  node = bm_render_shade (root, node, &face->shade,
                          face->width, face->height);

  return bm_render_warp (root, node, &matrix);
}

/**
 * bm_render_stack:
 * @root: the graph root
 * @background: a GEGL color string, or "none" for transparent
 * @layers: (array length=n_layers): layers, bottom first
 * @n_layers: number of layers
 * @width: canvas width
 * @height: canvas height
 *
 * Composites @layers over a solid background with chained `gegl:over`
 * nodes and crops to the canvas.
 *
 * Returns: (transfer none): the crop node
 */
static GeglNode *
bm_render_stack (GeglNode    *root,
                 const char  *background,
                 GeglNode   **layers,
                 int          n_layers,
                 int          width,
                 int          height)
{
  GeglColor *color;
  GeglNode *last;
  GeglNode *node;
  int i;

  if (g_strcmp0 (background, "none") == 0)
    color = bm_render_gray (0.0, 0.0);
  else
    color = gegl_color_new (background);

  node = gegl_node_new_child (root,
                              "operation", "gegl:color",
                              "value", color,
                              NULL);
  g_object_unref (color);
  last = gegl_node_new_child (root,
                              "operation", "gegl:crop",
                              "x", 0.0,
                              "y", 0.0,
                              "width", (double) width,
                              "height", (double) height,
                              NULL);
  gegl_node_link (node, last);

  for (i = 0; i < n_layers; i++)
    {
      node = gegl_node_new_child (root,
                                  "operation", "gegl:over",
                                  NULL);
      gegl_node_link (last, node);
      gegl_node_connect (node, "aux", layers[i], "output");
      last = node;
    }

  node = gegl_node_new_child (root,
                              "operation", "gegl:crop",
                              "x", 0.0,
                              "y", 0.0,
                              "width", (double) width,
                              "height", (double) height,
                              NULL);
  gegl_node_link (last, node);

  return node;
}

/**
 * bm_render_save:
 * @root: the graph root
 * @input: the final node
 * @path: output file; the format follows the extension
 *
 * Adds `gegl:save` and processes the graph.
 */
static void
bm_render_save (GeglNode   *root,
                GeglNode   *input,
                const char *path)
{
  GeglNode *save;

  save = gegl_node_new_child (root,
                              "operation", "gegl:save",
                              "path", path,
                              NULL);
  gegl_node_link (input, save);
  gegl_node_process (save);
}


/* Stratum 4: App
 *
 * Everything specific to books: options, command line, and the mapping
 * from cover, spine and page artwork onto the sides of a box. Replace
 * this stratum for other box-shaped products.
 */

typedef struct _BmOptions BmOptions;

struct _BmOptions
{
  char   *cover;
  char   *spine;
  char   *back;             /* optional; plain board when unset         */
  char   *output;

  double  angle;            /* degrees; 0 front, 90 spine, 180 back     */
  double  focal;            /* pixels; larger is flatter                */
  double  camera_height;    /* eye above the floor, fraction of height  */
  int     canvas_width;
  int     canvas_height;
  int     book_height;      /* rendered height of the hinge edge        */
  char   *background;       /* GEGL color string, or "none"             */
  double  light_ambient;    /* brightness of faces turned from light    */
  double  light_diffuse;    /* added brightness for faces facing it     */
  double  depth_falloff;    /* 0..1 brightness at the book's far depth  */
  double  page_color[3];    /* linear RGB of the page edges             */
  double  board_color[3];   /* linear RGB of a back cover without art   */
  double  page_lines;       /* 0..1 strength of streaked page lines     */
  double  page_grain;       /* 0..1 strength of fine paper grain        */
  char   *texture;          /* optional texture file for cover + spine  */
  double  texture_strength; /* 0..1                                     */
  double  cover_grain;      /* 0..1 noise grain on cover + spine        */
  int     seed;             /* base seed for generated noise            */
  double  crease;           /* 0..1 depth of the cover's hinge crease   */
  double  crease_offset;    /* distance from spine, fraction of height  */
  double  crease_width;     /* groove width, fraction of height         */

  double  shadow_direction; /* degrees from straight back, + is left    */
  double  shadow_length;    /* fraction of book height                  */
  double  shadow_opacity;
  double  shadow_blur;
  int     shadow_steps;
  double  contact_opacity;
  double  contact_blur;
};

/* String options stay %NULL here and get their defaults after parsing, so
 * `GOptionContext` never overwrites an allocated default.
 */
static void
bm_options_init_defaults (BmOptions *options)
{
  *options = (BmOptions) {
    .angle            = 25.0,
    .focal            = 6000.0,
    .camera_height    = 0.5,
    .canvas_width     = 2400,
    .canvas_height    = 2400,
    .book_height      = 1800,
    .light_ambient    = 0.70,
    .light_diffuse    = 0.60,
    .depth_falloff    = 0.82,
    .page_color       = { 0.86, 0.83, 0.76 },
    .board_color      = { 0.80, 0.80, 0.80 },
    .page_lines       = 0.18,
    .page_grain       = 0.05,
    .texture_strength = 0.15,
    .cover_grain      = 0.0,
    .seed             = 1,
    .crease           = 0.35,
    .crease_offset    = 0.03,
    .crease_width     = 0.006,
    .shadow_direction = 40.0,
    .shadow_length    = 0.55,
    .shadow_opacity   = 0.55,
    .shadow_blur      = 30.0,
    .shadow_steps     = 10,
    .contact_opacity  = 0.85,
    .contact_blur     = 4.0,
  };
}

static void
bm_options_clear (BmOptions *options)
{
  g_clear_pointer (&options->cover, g_free);
  g_clear_pointer (&options->spine, g_free);
  g_clear_pointer (&options->back, g_free);
  g_clear_pointer (&options->output, g_free);
  g_clear_pointer (&options->background, g_free);
  g_clear_pointer (&options->texture, g_free);
}

/**
 * bm_options_parse:
 * @options: a `BmOptions` with defaults already set
 * @argc: (inout): pointer to the argument count
 * @argv: (inout): pointer to the argument vector
 * @error: return location for a `GError`
 *
 * Parses the command line with `GOptionContext`.
 *
 * Returns: %TRUE on success
 */
static gboolean
bm_options_parse (BmOptions   *options,
                  int         *argc,
                  char      ***argv,
                  GError     **error)
{
  const GOptionEntry entries[] = {
    { "cover", 'c', 0, G_OPTION_ARG_FILENAME, &options->cover,
      "Front cover image", "FILE" },
    { "spine", 's', 0, G_OPTION_ARG_FILENAME, &options->spine,
      "Spine image", "FILE" },
    { "back", 'k', 0, G_OPTION_ARG_FILENAME, &options->back,
      "Back cover image (default: plain board)", "FILE" },
    { "output", 'o', 0, G_OPTION_ARG_FILENAME, &options->output,
      "Output file (default: mockup.png)", "FILE" },
    { "angle", 'a', 0, G_OPTION_ARG_DOUBLE, &options->angle,
      "Rotation in degrees: 0 front, 90 spine, 180 back, 155 back with "
      "spine (default: 25)", "DEG" },
    { "focal", 'f', 0, G_OPTION_ARG_DOUBLE, &options->focal,
      "Focal length; larger is flatter (default: 6000)", "PX" },
    { "width", 'W', 0, G_OPTION_ARG_INT, &options->canvas_width,
      "Canvas width (default: 2400)", "PX" },
    { "height", 'H', 0, G_OPTION_ARG_INT, &options->canvas_height,
      "Canvas height (default: 2400)", "PX" },
    { "book-height", 'b', 0, G_OPTION_ARG_INT, &options->book_height,
      "Book height on the canvas (default: 1800)", "PX" },
    { "camera-height", 'e', 0, G_OPTION_ARG_DOUBLE, &options->camera_height,
      "Eye height above the floor, as a fraction of book height; "
      "above 1 looks down on the top (default: 0.5)", "N" },
    { "background", 'g', 0, G_OPTION_ARG_STRING, &options->background,
      "Background color, or \"none\" (default: white)", "COLOR" },
    { "shadow-direction", 0, 0, G_OPTION_ARG_DOUBLE,
      &options->shadow_direction,
      "Shadow direction from straight back, + is left (default: 40)",
      "DEG" },
    { "shadow-length", 0, 0, G_OPTION_ARG_DOUBLE, &options->shadow_length,
      "Shadow length as a fraction of height (default: 0.55)", "N" },
    { "shadow-opacity", 0, 0, G_OPTION_ARG_DOUBLE, &options->shadow_opacity,
      "Cast shadow opacity (default: 0.55)", "N" },
    { "page-lines", 0, 0, G_OPTION_ARG_DOUBLE, &options->page_lines,
      "Strength of page-edge lines; 0 disables (default: 0.18)", "N" },
    { "page-grain", 0, 0, G_OPTION_ARG_DOUBLE, &options->page_grain,
      "Strength of paper grain on the pages; 0 disables (default: 0.05)",
      "N" },
    { "texture", 't', 0, G_OPTION_ARG_FILENAME, &options->texture,
      "Texture file applied to the cover and spine", "FILE" },
    { "texture-strength", 0, 0, G_OPTION_ARG_DOUBLE,
      &options->texture_strength,
      "Strength of --texture (default: 0.15)", "N" },
    { "cover-grain", 0, 0, G_OPTION_ARG_DOUBLE, &options->cover_grain,
      "Strength of noise grain on the cover and spine (default: 0)", "N" },
    { "seed", 0, 0, G_OPTION_ARG_INT, &options->seed,
      "Base seed for generated noise (default: 1)", "N" },
    { "crease", 0, 0, G_OPTION_ARG_DOUBLE, &options->crease,
      "Depth of the hinge creases on the covers; 0 disables "
      "(default: 0.35)", "N" },
    { "crease-offset", 0, 0, G_OPTION_ARG_DOUBLE, &options->crease_offset,
      "Crease distance from the spine, as a fraction of book height "
      "(default: 0.03)", "N" },
    { "crease-width", 0, 0, G_OPTION_ARG_DOUBLE, &options->crease_width,
      "Crease width, as a fraction of book height (default: 0.006)", "N" },
    { NULL }
  };
  g_autoptr (GOptionContext) context = NULL;

  context = g_option_context_new ("- render an upright book mockup");
  g_option_context_add_main_entries (context, entries, NULL);

  if (!g_option_context_parse (context, argc, argv, error))
    return FALSE;

  if (options->cover == NULL || options->spine == NULL)
    {
      g_set_error_literal (error, G_OPTION_ERROR, G_OPTION_ERROR_FAILED,
                           "both --cover and --spine are required");
      return FALSE;
    }

  if (options->shadow_steps > BM_MAX_BANDS)
    options->shadow_steps = BM_MAX_BANDS;

  if (options->output == NULL)
    options->output = g_strdup ("mockup.png");
  if (options->background == NULL)
    options->background = g_strdup ("white");

  return TRUE;
}

/**
 * bm_book_image_size:
 * @root: the graph root
 * @path: image file
 * @width: (out): return location for the width
 * @height: (out): return location for the height
 * @error: return location for a `GError`
 *
 * Reads an image's pixel size through `gegl:load`.
 *
 * Returns: %TRUE if the image could be read
 */
static gboolean
bm_book_image_size (GeglNode    *root,
                    const char  *path,
                    int         *width,
                    int         *height,
                    GError     **error)
{
  GeglNode *load;
  GeglRectangle box;

  /* gegl:load substitutes a placeholder for missing files, so check
   * first.
   */
  if (!g_file_test (path, G_FILE_TEST_IS_REGULAR))
    {
      g_set_error (error, G_FILE_ERROR, G_FILE_ERROR_NOENT,
                   "no such file “%s”", path);
      return FALSE;
    }

  load = gegl_node_new_child (root,
                              "operation", "gegl:load",
                              "path", path,
                              NULL);
  box = gegl_node_get_bounding_box (load);

  if (box.width <= 0 || box.height <= 0)
    {
      g_set_error (error, G_FILE_ERROR, G_FILE_ERROR_FAILED,
                   "could not read image “%s”", path);
      return FALSE;
    }

  *width = box.width;
  *height = box.height;

  return TRUE;
}

/**
 * bm_book_light:
 * @options: the render options
 *
 * Returns: the scene light described by @options
 */
static BmLight
bm_book_light (const BmOptions *options)
{
  return (BmLight) {
    .azimuth = options->shadow_direction,
    .length  = options->shadow_length,
    .ambient = options->light_ambient,
    .diffuse = options->light_diffuse,
  };
}

/**
 * bm_book_add_page_textures:
 * @face: a page-edge face
 * @options: the render options
 * @angle: direction of the page lines in the face's source space
 * @length: streak length, in source units
 * @seed: seed for the first noise layer
 *
 * Adds streaked page lines and fine paper grain to @face.
 */
static void
bm_book_add_page_textures (BmFace          *face,
                           const BmOptions *options,
                           double           angle,
                           double           length,
                           guint            seed)
{
  if (options->page_lines > 0.0)
    bm_face_add_texture (face, &(BmTexture) {
                           .kind     = BM_TEXTURE_NOISE,
                           .strength = options->page_lines,
                           .streak   = length,
                           .angle    = angle,
                           .seed     = seed,
                         });
  if (options->page_grain > 0.0)
    bm_face_add_texture (face, &(BmTexture) {
                           .kind     = BM_TEXTURE_NOISE,
                           .strength = options->page_grain,
                           .seed     = seed + 1,
                         });
}

/**
 * bm_book_build_scene:
 * @scene: (out caller-allocates): the scene to fill in
 * @root: the graph root, used to read the artwork sizes
 * @options: the render options
 * @error: return location for a `GError`
 *
 * Builds the book as a six-sided box: the front cover, the back cover
 * (or plain board), the spine on the left side, and page edges on the
 * right (fore-edge), top and bottom. The box is sized from the front
 * cover's aspect ratio and the spine's width at the book height, then
 * rotated, lit and centered.
 *
 * Returns: %TRUE on success
 */
static gboolean
bm_book_build_scene (BmScene          *scene,
                     GeglNode         *root,
                     const BmOptions  *options,
                     GError          **error)
{
  double h = options->book_height;
  int cover_w, cover_h, spine_w, spine_h;
  double w, d;
  BmFace front, back, spine, fore_edge, top, bottom;
  const BmFace *sides[BM_BOX_N_SIDES];
  BmLight light;

  if (!bm_book_image_size (root, options->cover, &cover_w, &cover_h, error) ||
      !bm_book_image_size (root, options->spine, &spine_w, &spine_h, error))
    return FALSE;

  if (options->back != NULL)
    {
      int bw, bh;

      if (!bm_book_image_size (root, options->back, &bw, &bh, error))
        return FALSE;
    }

  /* Both pieces are scaled to the book height, keeping their aspect
   * ratios, so spine thickness stays true to the trim size as long as both
   * came from the same cover wrap at the same resolution. The back cover
   * is stretched to the front's size.
   */
  w = cover_w * h / cover_h;
  d = spine_w * h / spine_h;

  front = (BmFace) { .image = options->cover, .shade = { 1.0, 1.0 } };
  spine = (BmFace) { .image = options->spine, .shade = { 1.0, 1.0 } };
  back = (BmFace) {
    .image = options->back,
    .color = { options->board_color[0], options->board_color[1],
               options->board_color[2] },
    .shade = { 1.0, 1.0 },
  };
  fore_edge = (BmFace) {
    .color = { options->page_color[0], options->page_color[1],
               options->page_color[2] },
    .shade = { 1.0, 1.0 },
  };
  top = fore_edge;
  bottom = fore_edge;

  /* Page edges: each sheet's edge is a line along the face, stacked
   * through the book's thickness. On the top and bottom the lines run
   * along the source x axis (parallel to the covers); on the fore-edge
   * they run along the source y axis.
   */
  bm_book_add_page_textures (&top, options, 0.0, w / 4.0, options->seed);
  bm_book_add_page_textures (&bottom, options, 0.0, w / 4.0,
                             options->seed + 2);
  bm_book_add_page_textures (&fore_edge, options, 90.0, h / 4.0,
                             options->seed + 4);

  /* Cover stock: the same material on both covers and the spine */
  if (options->texture != NULL)
    {
      int tw, th;
      BmTexture stock = {
        .kind     = BM_TEXTURE_IMAGE,
        .image    = options->texture,
        .strength = options->texture_strength,
      };

      if (!bm_book_image_size (root, options->texture, &tw, &th, error))
        return FALSE;

      bm_face_add_texture (&front, &stock);
      bm_face_add_texture (&back, &stock);
      bm_face_add_texture (&spine, &stock);
    }
  if (options->cover_grain > 0.0)
    {
      BmTexture grain = {
        .kind     = BM_TEXTURE_NOISE,
        .strength = options->cover_grain,
        .seed     = options->seed + 6,
      };

      bm_face_add_texture (&front, &grain);
      grain.seed++;
      bm_face_add_texture (&back, &grain);
      grain.seed++;
      bm_face_add_texture (&spine, &grain);
    }

  /* Hinge creases: score lines running parallel to the spine a short
   * distance in from it. The spine is on the front's left and, seen from
   * behind, on the back's right.
   */
  if (options->crease > 0.0)
    {
      BmTexture crease = {
        .kind     = BM_TEXTURE_CREASE,
        .strength = options->crease,
        .offset   = options->crease_offset * h,
        .width    = options->crease_width * h,
      };

      bm_face_add_texture (&front, &crease);
      crease.offset = w - crease.offset;
      bm_face_add_texture (&back, &crease);
    }

  sides[BM_BOX_FRONT] = &front;
  sides[BM_BOX_BACK] = &back;
  sides[BM_BOX_LEFT] = &spine;
  sides[BM_BOX_RIGHT] = &fore_edge;
  sides[BM_BOX_TOP] = &top;
  sides[BM_BOX_BOTTOM] = &bottom;

  light = bm_book_light (options);

  *scene = (BmScene) { .camera = { .focal = options->focal } };
  bm_scene_add_box (scene, sides, w, d, h, options->angle,
                    options->camera_height * h);
  bm_scene_light (scene, &light, options->depth_falloff);
  bm_scene_center (scene, options->canvas_width, options->canvas_height);

  return TRUE;
}

int
main (int    argc,
      char **argv)
{
  g_autoptr (GError) error = NULL;
  g_autoptr (GeglNode) root = NULL;
  BmOptions options;
  BmScene scene;
  BmLight light;
  BmBand cast_bands[BM_MAX_BANDS];
  BmBand contact_band;
  GeglNode *masks[2];
  GeglNode *layers[1 + BM_MAX_FACES];
  GeglNode *result;
  int n_cast;
  int n_layers = 0;
  int status = 0;
  int i;

  bm_options_init_defaults (&options);

  if (!bm_options_parse (&options, &argc, &argv, &error))
    {
      g_printerr ("book-mockup: %s\n", error->message);
      bm_options_clear (&options);
      return 1;
    }

  gegl_init (NULL, NULL);
  root = gegl_node_new ();

  if (!bm_book_build_scene (&scene, root, &options, &error))
    {
      g_printerr ("book-mockup: %s\n", error->message);
      status = 1;
      goto out;
    }

  /* Floor shadow: a long, soft cast shadow plus a tight contact line */
  light = bm_book_light (&options);
  n_cast = bm_scene_shadow_bands (&scene, &light, 1.0, options.shadow_steps,
                                  options.shadow_opacity, cast_bands);
  bm_scene_shadow_bands (&scene, &light, 0.04, 1,
                         options.contact_opacity, &contact_band);

  masks[0] = bm_render_mask (root, cast_bands, n_cast, options.shadow_blur,
                             options.canvas_width, options.canvas_height);
  masks[1] = bm_render_mask (root, &contact_band, 1, options.contact_blur,
                             options.canvas_width, options.canvas_height);
  layers[n_layers++] = bm_render_shadow (root, masks, G_N_ELEMENTS (masks),
                                         options.canvas_width,
                                         options.canvas_height);

  /* Faces that point toward the camera. The box is convex, so visible
   * faces never overlap and their order does not matter.
   */
  for (i = 0; i < scene.n_faces; i++)
    {
      GeglNode *face;

      if (!bm_face_is_visible (&scene.faces[i], &scene.camera))
        continue;

      face = bm_render_face (root, &scene.faces[i], &scene.camera);
      if (face != NULL)
        layers[n_layers++] = face;
    }

  result = bm_render_stack (root, options.background, layers, n_layers,
                            options.canvas_width, options.canvas_height);
  bm_render_save (root, result, options.output);

out:
  /* The graph must be released before gegl_exit(); clearing it here
   * leaves the autoptr with nothing to do.
   */
  g_clear_object (&root);
  gegl_exit ();
  bm_options_clear (&options);

  return status;
}
