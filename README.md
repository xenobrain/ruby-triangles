# ruby-triangles

A single-file Ruby library for 2D computational geometry. Provides Delaunay triangulation, convex hull, polygon triangulation with holes, visibility polygons, constrained Delaunay triangulation, and polygon decomposition/shattering utilities.

All algorithms use **flat coordinate arrays** — coordinates are stored as `[x1, y1, x2, y2, ...]` rather than nested point objects, which keeps memory overhead low and makes it easy to pass data between the different algorithms.

## Algorithms

| Algorithm | Method | Notes |
|-----------|--------|-------|
| Delaunay triangulation | `Triangles.delaunay` | Randomized incremental, O(n log n) |
| Convex hull | `Triangles.convex_hull` | QuickHull, O(n log n) |
| Polygon triangulation | `Triangles.earcut` | Ear-clipping with holes, O(n²) |
| Convex decomposition | `Triangles.decompose_to_convex_hulls` | Bayazit algorithm |
| Constrained Delaunay | `Triangles::TriVis.cdt_init` | Edge-constrained CDT |
| Visibility polygon | `Triangles::TriVis.triangular_expansion` | Triangular expansion |
| Voronoi shattering | `Triangles.shatter_polygon` | Worley-noise based |
| Robust predicates | `Triangles::Geometry` | Shewchuk adaptive arithmetic |

## Installation

Copy `triangles.rb` into your project and require it:

```ruby
require_relative 'triangles'
```

No gems required. The file is self-contained and has no external dependencies.

## Usage

### Delaunay Triangulation

```ruby
coords = [0, 0, 100, 0, 50, 87, 100, 100, 0, 100]
result = Triangles.delaunay(coords)

result[:triangles]   # => [i, j, k, ...] — vertex indices, one triple per triangle
result[:half_edges]  # => half-edge adjacency array; -1 means hull edge
result[:hull]        # => convex hull vertex indices in CCW order
result[:coords]      # => the original flat coordinate array
```

Collinear or fewer than three points return an empty triangulation.

### Convex Hull

```ruby
coords = [0, 0, 100, 0, 50, 50, 100, 100, 0, 100]
hull = Triangles.convex_hull(coords)
# => [0, 1, 3, 4] — vertex indices in CCW order

hull = Triangles.convex_hull(coords, tol = 0.5)  # optional collinearity tolerance
```

### Polygon Triangulation (Earcut)

Triangulates simple polygons and polygons with holes. Coordinates are flat; holes are specified by their starting vertex index.

```ruby
# Simple polygon
verts   = [0, 0, 10, 0, 10, 10, 0, 10]
indices = Triangles.earcut(verts)
# => [1, 0, 3, 3, 2, 1]

# Polygon with a hole
outer  = [0, 0, 100, 0, 100, 100, 0, 100]
hole   = [20, 20, 20, 80, 80, 80, 80, 20]
verts  = outer + hole
holes  = [outer.length / 2]   # hole starts at vertex index 4
indices = Triangles.earcut(verts, holes)

# 3D coordinates (only x, y are used for triangulation)
indices = Triangles.earcut(verts, nil, 3)
```

Returns a flat array of vertex indices: every three consecutive values form one triangle.

### Polygon Utilities

All utility methods use the flat `[x1, y1, x2, y2, ...]` coordinate convention.

```ruby
coords = [0, 0, 10, 0, 10, 10, 0, 10]

Triangles.polygon_area(coords)         # => -100.0 (negative = CW, positive = CCW)
Triangles.centroid(coords)             # => [5.0, 5.0]
Triangles.bounding_box(coords)         # => [min_x, min_y, max_x, max_y]
Triangles.convex?(coords)              # => true
Triangles.point_in_poly?(5, 5, coords) # => true

Triangles.ensure_clockwise(coords)     # reverses if CCW
Triangles.ensure_ccw(coords)           # reverses if CW
```

### Convex Decomposition (Bayazit)

Decomposes a concave polygon into convex sub-polygons. All output parts are guaranteed convex and preserve the input winding order.

```ruby
# L-shape (concave)
l_shape = [0, 0, 10, 0, 10, 5, 5, 5, 5, 10, 0, 10]
parts   = Triangles.decompose_to_convex_hulls(l_shape)
# => [[0, 0, 10, 0, 10, 5, ...], [5, 5, 5, 10, 0, 10, ...]]

parts.each { |p| puts Triangles.convex?(p) }  # => true for all parts
```

Convex input is returned unchanged in a single-element array.

### Constrained Delaunay Triangulation (CDT)

Adds hard edge constraints to a Delaunay triangulation, ensuring specified edges appear in the output.

```ruby
coords = [150, 50, 50, 200, 150, 350, 250, 200].flatten
del = Triangles.delaunay(coords)
del[:coords] = coords

# Constrain the diagonal edge between vertices 0 and 2
con = Triangles::TriVis.cdt_init(del, [[0, 2]])

# Query whether a specific half-edge is constrained
Triangles::TriVis.cdt_edge_constrained?(con, edge_index)  # => true/false

# Add individual constraints after init
Triangles::TriVis.cdt_constrain_edge(con, p1_index, p2_index)
```

### Visibility Polygon

Computes the polygon of all points visible from a query position inside a triangulated scene. Edge constraints act as opaque walls.

```ruby
# Build scene: triangulate points and constrain boundary edges
points = [[53, 98], [5, 201], [194, 288], [392, 148], [413, 43]]
coords = points.flatten
del = Triangles.delaunay(coords)
del[:coords] = coords

con = Triangles::TriVis.cdt_init(del)
n = coords.length / 2
n.times { |i| Triangles::TriVis.cdt_constrain_edge(con, i, (i + 1) % n) }

# Query: what is visible from (200, 150)?
obstructs = ->(edg) { Triangles::TriVis.cdt_edge_constrained?(con, edg) }
segments  = Triangles::TriVis.triangular_expansion(del, 200, 150, obstructs)
# => [[x1, y1, x2, y2], ...] — CCW-ordered visibility boundary segments

# Restrict to a viewing cone (optional)
segments = Triangles::TriVis.triangular_expansion(
  del, qx, qy, obstructs,
  left_x, left_y, right_x, right_y
)
```

`triangular_expansion` returns an empty array when the query point lies exactly on a vertex or outside the triangulation.

### Polygon Shattering

Breaks a polygon into Voronoi-like shards, useful for destruction effects.

```ruby
polygon = [0, 0, 200, 0, 200, 200, 0, 200]
shards  = Triangles.shatter_polygon(
  polygon,
  focus:        [100, 100],  # impact point (gets denser sites)
  cell_size:    40,           # grid spacing for Voronoi sites
  jitter:       0.4,          # randomness in site placement (0–1)
  min_area:     20.0,         # discard shards smaller than this
  focus_radius: 60,           # radius of dense zone around focus
  focus_extra_sites: 8,       # extra random sites near focus
  seed:         42            # for reproducibility
)
# => [[x1, y1, x2, y2, ...], ...] — one flat polygon per shard
```

## Reference

### `Triangles::Geometry` — Robust Predicates

```ruby
Triangles::Geometry.triangle_signed_area(ax, ay, bx, by, cx, cy)
# => positive (CCW), negative (CW), or ~0 (collinear)

Triangles::Geometry.robust_orient2d(ax, ay, bx, by, cx, cy)
# => exact sign using Shewchuk adaptive arithmetic

Triangles::Geometry.robust_incircle(ax, ay, bx, by, cx, cy, dx, dy)
# => > 0 if d inside circumcircle of CCW triangle (abc), <= 0 otherwise

Triangles::Geometry.segments_intersect?(p1x, p1y, p2x, p2y, p3x, p3y, p4x, p4y)
# => true if segments [p1,p2] and [p3,p4] properly intersect

Triangles::Geometry.point_on_segment?(ax, ay, bx, by, px, py)
# => true if p lies on segment [a,b] (bounding-box check)

Triangles::Geometry.dist_squared(ax, ay, bx, by)
# => squared Euclidean distance (avoids sqrt)
```

### `Triangles::Shatter` — Low-level Polygon Operations

```ruby
# Sutherland–Hodgman half-plane clipping
# Keeps the side where dot(point, [nx,ny]) <= threshold
Triangles::Shatter.clip_polygon_half_plane(coords, nx, ny, threshold)

# Test whether vertex i of a flat polygon is a reflex (concave) vertex
Triangles::Shatter.bayazit_reflex?(i, coords)

# Line–line intersection point (returns [0, 0] for parallel lines)
Triangles::Shatter.line_intersect(p1x, p1y, p2x, p2y, p3x, p3y, p4x, p4y)

# Generate Worley noise sites inside a polygon for shattering
Triangles::Shatter.generate_worley_sites(coords, cell_size:, jitter:,
  focus: nil, seed: nil, focus_radius: nil, focus_extra_sites: 3)
# => { inside: [[x,y],...], all: [[x,y],...], inside_flags: [true/false,...] }
```

## Running Tests

```bash
ruby tests/test_all.rb
```

Individual test suites:

```bash
ruby tests/test_triangles.rb      # Delaunay triangulation
ruby tests/test_quickhull.rb      # Convex hull
ruby tests/test_earcut.rb         # Polygon triangulation
ruby tests/test_bayazit.rb        # Convex decomposition
ruby tests/test_shatter_helpers.rb # Utility functions
ruby tests/test_constrain.rb      # Constrained Delaunay
ruby tests/test_trivis.rb         # Visibility polygons
```

Tests require the `minitest` gem (included in Ruby's standard library).

## Design Notes

- **Single file**: The entire library lives in `triangles.rb`. No gem, no bundler required.
- **Flat arrays**: All coordinate data uses `[x1, y1, x2, y2, ...]`. This avoids object allocation overhead and makes the format easy to share between algorithms.
- **Numerically robust**: Orientation and in-circle tests use Shewchuk's adaptive exact arithmetic, eliminating the "wobbles" that plague naive floating-point geometry.
- **No mutation of input**: Methods return new arrays; input coordinates are not modified.
- **Public domain**: Licensed under the Unlicense — do whatever you want with it.

## License

This is free and unencumbered software released into the public domain. See [LICENSE](LICENSE) for details.
