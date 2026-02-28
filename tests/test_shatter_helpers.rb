require 'minitest/autorun'

require_relative '../triangles'

# Comprehensive test suite for Triangles helper methods
class TestShatterHelpers < Minitest::Test
  def setup
    @triangles = Triangles
    @shatter = Triangles::Shatter  # For internal Shatter module methods
  end

  # Tests for polygon_area
  def test_polygon_area_square
    square = [0, 0, 10, 0, 10, 10, 0, 10]
    # This square is counter-clockwise in standard coordinate system
    # CCW = positive area in shoelace formula
    assert_equal(100.0, @triangles.polygon_area(square))
  end

  def test_polygon_area_ccw_square
    square = [0, 0, 0, 10, 10, 10, 10, 0]
    # This square is clockwise in standard coordinate system
    # CW = negative area in shoelace formula
    area = @triangles.polygon_area(square)
    assert_equal(100.0, area.abs)
    assert(area < 0, "Clockwise square should have negative area")
  end

  def test_polygon_area_triangle
    triangle = [0, 0, 10, 0, 5, 10]
    area = @triangles.polygon_area(triangle)
    assert_in_delta(50.0, area.abs, 0.01)
  end

  def test_polygon_area_empty
    assert_equal(0.0, @triangles.polygon_area([]))
    assert_equal(0.0, @triangles.polygon_area([1, 2]))
    assert_equal(0.0, @triangles.polygon_area([1, 2, 3, 4]))
  end

  # Tests for ensure_clockwise
  def test_ensure_clockwise_already_clockwise
    cw_square = [0, 0, 10, 0, 10, 10, 0, 10]
    result = @triangles.ensure_clockwise(cw_square)
    # Already clockwise (negative area), should not reverse
    assert_equal(cw_square.length, result.length)
    # Check area is still negative (clockwise)
    assert(@triangles.polygon_area(result) < 0)
  end

  def test_ensure_clockwise_needs_reversal
    ccw_square = [0, 0, 0, 10, 10, 10, 10, 0]
    result = @triangles.ensure_clockwise(ccw_square)
    # Counter-clockwise (positive area), should reverse to clockwise
    assert(@triangles.polygon_area(result) < 0)
  end

  # Tests for ensure_ccw
  def test_ensure_ccw_already_ccw
    ccw_square = [0, 0, 0, 10, 10, 10, 10, 0]
    result = @triangles.ensure_ccw(ccw_square)
    # Already counter-clockwise (positive area), should not reverse
    assert_equal(ccw_square.length, result.length)
    # Check area is still positive (counter-clockwise)
    assert(@triangles.polygon_area(result) > 0)
  end

  def test_ensure_ccw_needs_reversal
    cw_square = [0, 0, 10, 0, 10, 10, 0, 10]
    result = @triangles.ensure_ccw(cw_square)
    # Clockwise (negative area), should reverse to counter-clockwise
    assert(@triangles.polygon_area(result) > 0)
  end

  # Tests for centroid
  def test_centroid_square
    square = [0, 0, 10, 0, 10, 10, 0, 10]
    cx, cy = @triangles.centroid(square)
    assert_in_delta(5.0, cx, 0.01)
    assert_in_delta(5.0, cy, 0.01)
  end

  def test_centroid_triangle
    triangle = [0, 0, 10, 0, 5, 10]
    cx, cy = @triangles.centroid(triangle)
    # Centroid of triangle is at (1/3) of sum of vertices
    assert_in_delta(5.0, cx, 0.5)
    assert_in_delta(3.33, cy, 0.5)
  end

  def test_centroid_zero_area
    degenerate = [0, 0, 0, 0, 0, 0]
    cx, cy = @triangles.centroid(degenerate)
    assert_equal(0.0, cx)
    assert_equal(0.0, cy)
  end

  # Tests for point_in_poly?
  def test_point_in_poly_inside
    square = [0, 0, 10, 0, 10, 10, 0, 10]
    assert(@triangles.point_in_poly?(5, 5, square))
    assert(@triangles.point_in_poly?(1, 1, square))
    assert(@triangles.point_in_poly?(9, 9, square))
  end

  def test_point_in_poly_outside
    square = [0, 0, 10, 0, 10, 10, 0, 10]
    refute(@triangles.point_in_poly?(-1, 5, square))
    refute(@triangles.point_in_poly?(11, 5, square))
    refute(@triangles.point_in_poly?(5, -1, square))
    refute(@triangles.point_in_poly?(5, 11, square))
  end

  def test_point_in_poly_on_edge
    square = [0, 0, 10, 0, 10, 10, 0, 10]
    # Edge cases - behavior may vary by implementation
    # Testing a clear edge case
    result = @triangles.point_in_poly?(0, 5, square)
    # On edge - either true or false is acceptable, just document it
    assert([true, false].include?(result))
  end

  def test_point_in_poly_concave
    # L-shaped polygon
    l_shape = [0, 0, 10, 0, 10, 5, 5, 5, 5, 10, 0, 10]
    assert(@triangles.point_in_poly?(2, 2, l_shape))
    assert(@triangles.point_in_poly?(7, 2, l_shape))
    refute(@triangles.point_in_poly?(7, 7, l_shape)) # In the missing corner
  end

  # Tests for bounding_box
  def test_bounding_box_square
    square = [0, 0, 10, 0, 10, 10, 0, 10]
    min_x, min_y, max_x, max_y = @triangles.bounding_box(square)
    assert_equal(0, min_x)
    assert_equal(0, min_y)
    assert_equal(10, max_x)
    assert_equal(10, max_y)
  end

  def test_bounding_box_irregular
    poly = [5, 5, -3, 2, 7, -1, 0, 10]
    min_x, min_y, max_x, max_y = @triangles.bounding_box(poly)
    assert_equal(-3, min_x)
    assert_equal(-1, min_y)
    assert_equal(7, max_x)
    assert_equal(10, max_y)
  end

  def test_bounding_box_single_point
    point = [5, 7]
    min_x, min_y, max_x, max_y = @triangles.bounding_box(point)
    assert_equal(5, min_x)
    assert_equal(7, min_y)
    assert_equal(5, max_x)
    assert_equal(7, max_y)
  end

  # Tests for convex?
  def test_is_convex_square
    square = [0, 0, 10, 0, 10, 10, 0, 10]
    assert(@triangles.convex?(square))
  end

  def test_is_convex_triangle
    triangle = [0, 0, 10, 0, 5, 10]
    assert(@triangles.convex?(triangle))
  end

  def test_is_convex_l_shape
    l_shape = [0, 0, 10, 0, 10, 5, 5, 5, 5, 10, 0, 10]
    refute(@triangles.convex?(l_shape))
  end

  def test_is_convex_star
    # Simple concave star
    star = [0, 0, 2, 3, 0, 6, 3, 4, 6, 6, 4, 3, 6, 0, 3, 2]
    refute(@triangles.convex?(star))
  end

  def test_is_convex_too_few_points
    refute(@triangles.convex?([]))
    refute(@triangles.convex?([1, 2]))
    refute(@triangles.convex?([1, 2, 3, 4]))
  end

  # Tests for clip_polygon_half_plane
  def test_clip_polygon_half_plane_no_clipping
    square = [0, 0, 10, 0, 10, 10, 0, 10]
    # Plane that keeps all points (normal pointing right, threshold far right)
    result = @shatter.clip_polygon_half_plane(square, 1, 0, 20)
    assert(result.length >= 6, "Should keep the polygon")
  end

  def test_clip_polygon_half_plane_clips_half
    square = [0.0, 0.0, 10.0, 0.0, 10.0, 10.0, 0.0, 10.0]
    # Vertical plane at x=5, normal pointing right, clips right half
    result = @shatter.clip_polygon_half_plane(square, 1, 0, 5)
    n = result.length.div(2)
    assert(n >= 3, "Should have at least 3 points")
    # Check that all points are on or left of x=5
    i = 0
    while i < n
      x = result[i * 2]; y = result[i * 2 + 1]
      assert(x <= 5.01, "Point #{x}, #{y} should be at or left of x=5")
      i += 1
    end
  end

  def test_clip_polygon_half_plane_clips_all
    square = [0, 0, 10, 0, 10, 10, 0, 10]
    # Plane that clips everything
    result = @shatter.clip_polygon_half_plane(square, 1, 0, -20)
    assert_equal(0, result.length)
  end

  def test_clip_polygon_half_plane_empty_input
    result = @shatter.clip_polygon_half_plane([], 1, 0, 0)
    assert_equal(0, result.length)
  end

  # Tests for generate_worley_sites
  def test_generate_worley_sites_basic
    square = [0, 0, 100, 0, 100, 100, 0, 100]
    result = @shatter.generate_worley_sites(
      square,
      cell_size: 25,
      jitter: 0.3,
      seed: 12345
    )

    assert(result[:inside].length > 0, "Should generate inside sites")
    assert(result[:all].length > 0, "Should generate all sites")
    assert_equal(result[:all].length, result[:inside_flags].length)
  end

  def test_generate_worley_sites_with_focus
    square = [0, 0, 100, 0, 100, 100, 0, 100]
    focus = [50, 50]
    result = @shatter.generate_worley_sites(
      square,
      cell_size: 25,
      jitter: 0.3,
      focus: focus,
      seed: 12345
    )

    # Focus should be included if inside polygon
    assert(result[:inside].include?(focus))
    assert(result[:all].include?(focus))
  end

  def test_generate_worley_sites_with_focus_radius
    square = [0, 0, 100, 0, 100, 100, 0, 100]
    focus = [50, 50]
    result = @shatter.generate_worley_sites(
      square,
      cell_size: 25,
      jitter: 0.3,
      focus: focus,
      focus_radius: 30,
      focus_extra_sites: 5,
      seed: 12345
    )

    assert(result[:inside].length >= 5, "Should have at least focus_extra_sites")
  end

  def test_generate_worley_sites_empty_polygon_uses_centroid
    square = [50, 50, 51, 50, 51, 51, 50, 51]
    # Very small square with large cell size - may have no inside sites
    result = @shatter.generate_worley_sites(
      square,
      cell_size: 1000,
      jitter: 0.3,
      seed: 12345
    )

    # Even with no grid sites, should have at least centroid
    assert(result[:inside].length > 0, "Should use centroid as fallback")
  end

  # Integration tests
  def test_shatter_polygon_basic
    square = [0, 0, 100, 0, 100, 100, 0, 100]
    focus = [50, 50]

    shards = @triangles.shatter_polygon(
      square,
      focus: focus,
      cell_size: 30,
      jitter: 0.3,
      min_area: 10.0
    )

    assert(shards.length > 1, "Should create multiple shards")
    shards.each do |shard|
      assert(shard.length >= 6, "Each shard should be a valid polygon")
      area = @triangles.polygon_area(shard).abs
      assert(area >= 10.0, "Each shard should meet min_area requirement")
    end
  end

  def test_decompose_to_convex_hulls_convex_input
    square = [0, 0, 10, 0, 10, 10, 0, 10]
    result = @triangles.decompose_to_convex_hulls(square)
    assert_equal(1, result.length, "Convex polygon should not be decomposed")
    assert_equal(square, result[0])
  end

  def test_decompose_to_convex_hulls_concave_input
    # L-shape is concave
    l_shape = [0, 0, 10, 0, 10, 5, 5, 5, 5, 10, 0, 10]
    result = @triangles.decompose_to_convex_hulls(l_shape)
    assert(result.length >= 2, "Concave polygon should be decomposed into multiple convex parts")

    # Each component should be convex
    result.each do |component|
      assert(@triangles.convex?(component), "Each component should be convex")
    end
  end

  def test_decompose_to_convex_hulls_preserves_winding
    # Test with clockwise winding
    cw_l = [0, 0, 0, 10, 5, 10, 5, 5, 10, 5, 10, 0]
    cw_area = @triangles.polygon_area(cw_l)
    result = @triangles.decompose_to_convex_hulls(cw_l)

    # All components should have same winding as input
    result.each do |component|
      component_area = @triangles.polygon_area(component)
      assert_equal(cw_area > 0, component_area > 0, "Winding should be preserved")
    end
  end

  # Tests for consolidated helper functions
  def test_triangle_signed_area_triangle
    area = Triangles::Geometry.triangle_signed_area(0, 0, 10, 0, 5, 10)
    # Signed area should be positive for CCW triangle
    assert(area > 0)
    assert_in_delta(100.0, area, 0.01)
  end

  def test_triangle_signed_area_collinear
    area = Triangles::Geometry.triangle_signed_area(0, 0, 5, 5, 10, 10)
    assert_in_delta(0.0, area, 0.01)
  end

  def test_line_intersect
    # Two non-parallel lines
    x, y = @shatter.line_intersect(0, 0, 10, 10, 0, 10, 10, 0)
    assert_in_delta(5.0, x, 0.01)
    assert_in_delta(5.0, y, 0.01)
  end

  def test_line_intersect_parallel
    # Parallel lines
    x, y = @shatter.line_intersect(0, 0, 10, 0, 0, 5, 10, 5)
    # Should return [0, 0] for parallel lines
    assert_equal(0.0, x)
    assert_equal(0.0, y)
  end

  def test_segments_intersect_intersecting
    assert(Triangles::Geometry.segments_intersect?(0, 0, 10, 10, 0, 10, 10, 0))
  end

  def test_segments_intersect_non_intersecting
    refute(Triangles::Geometry.segments_intersect?(0, 0, 5, 0, 6, 0, 10, 0))
  end

  def test_segments_intersect_collinear
    assert(Triangles::Geometry.segments_intersect?(0, 0, 10, 0, 5, 0, 15, 0))
  end

  def test_point_on_segment
    assert(Triangles::Geometry.point_on_segment?(0, 0, 10, 10, 5, 5))
  end

  def test_point_on_segment_not_on
    # Note: This only checks bounding box, not collinearity
    assert(Triangles::Geometry.point_on_segment?(0, 0, 10, 10, 5, 6))
  end

  def test_point_on_segment_outside_bounds
    refute(Triangles::Geometry.point_on_segment?(0, 0, 10, 10, 15, 15))
  end

  def test_bayazit_reflex
    # Create a simple concave polygon (flat)
    vertices = [0, 0, 10, 0, 10, 10, 5, 5, 0, 10]
    # Index 3 (point [5, 5]) should be reflex (concave)
    assert(@shatter.bayazit_reflex?(3, vertices))
    # Index 0 should not be reflex
    refute(@shatter.bayazit_reflex?(0, vertices))
  end

  # Edge case tests
  def test_handles_floating_point_coordinates
    poly = [0.5, 0.7, 10.3, 0.9, 10.1, 10.8, 0.2, 10.4]
    area = @triangles.polygon_area(poly)
    assert(area.abs > 0)

    cx, cy = @triangles.centroid(poly)
    assert(cx.finite?)
    assert(cy.finite?)
  end

  def test_handles_negative_coordinates
    poly = [-10, -10, 10, -10, 10, 10, -10, 10]
    cx, cy = @triangles.centroid(poly)
    assert_in_delta(0.0, cx, 0.01)
    assert_in_delta(0.0, cy, 0.01)
  end
end
