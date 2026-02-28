require 'minitest/autorun'

require_relative '../triangles'

class TestQuickHull < Minitest::Test
  def test_empty_input
    coords = []
    hull = Triangles.convex_hull(coords)
    assert_equal [], hull
  end

  def test_single_point
    coords = [5.0, 3.0]
    hull = Triangles.convex_hull(coords)
    assert_equal [0], hull
  end

  def test_two_points
    coords = [0.0, 0.0, 5.0, 5.0]
    hull = Triangles.convex_hull(coords)
    assert_equal [0, 1], hull
  end

  def test_three_points_triangle
    coords = [0.0, 0.0, 5.0, 0.0, 2.5, 5.0]
    hull = Triangles.convex_hull(coords)
    assert_equal 3, hull.length
    # All three points should be in the hull
    assert_equal [0, 1, 2].sort, hull.sort
  end

  def test_collinear_points
    # Points on a line
    coords = [0.0, 0.0, 1.0, 1.0, 2.0, 2.0, 3.0, 3.0]
    hull = Triangles.convex_hull(coords)
    # Should only include the endpoints
    assert_equal 2, hull.length
    assert_includes hull, 0  # leftmost point
    assert_includes hull, 3  # rightmost point
  end

  def test_square
    coords = [0.0, 0.0, 10.0, 0.0, 10.0, 10.0, 0.0, 10.0]
    hull = Triangles.convex_hull(coords)
    assert_equal 4, hull.length
    # All four corners should be in the hull
    assert_equal [0, 1, 2, 3].sort, hull.sort
  end

  def test_square_with_interior_point
    # Square with a point in the middle
    coords = [0.0, 0.0, 10.0, 0.0, 10.0, 10.0, 0.0, 10.0, 5.0, 5.0]
    hull = Triangles.convex_hull(coords)
    # Should only include the 4 corners, not the interior point
    assert_equal 4, hull.length
    assert_equal [0, 1, 2, 3].sort, hull.sort
  end

  def test_random_points
    # Generate random points
    srand(42)
    points = Array.new(20) { [rand * 100, rand * 100] }
    coords = points.flatten

    hull = Triangles.convex_hull(coords)

    # Hull should have at least 3 points for random data
    assert hull.length >= 3, "Hull should have at least 3 points"
    assert hull.length <= 20, "Hull should not have more points than input"

    # Verify all hull indices are valid
    hull.each do |idx|
      assert idx >= 0 && idx < 20, "Invalid hull index #{idx}"
    end

    # Verify hull is counter-clockwise by checking signed area
    area = 0.0
    hull.length.times do |i|
      x1, y1 = points[hull[i]]
      x2, y2 = points[hull[(i + 1) % hull.length]]
      area += x1 * y2 - x2 * y1
    end
    assert area > 0, "Hull should be counter-clockwise (positive area)"
  end

  def test_hull_order_is_correct
    # Pentagon-like shape
    coords = [0.0, 5.0, 3.0, 9.0, 8.0, 7.0, 7.0, 2.0, 2.0, 0.0]
    hull = Triangles.convex_hull(coords)

    # Verify points are in counter-clockwise order
    # Calculate signed area
    points = []
    i = 0
    while i < coords.length
      points << [coords[i], coords[i + 1]]
      i += 2
    end

    area = 0.0
    hull.length.times do |i|
      x1, y1 = points[hull[i]]
      x2, y2 = points[hull[(i + 1) % hull.length]]
      area += x1 * y2 - x2 * y1
    end

    assert area > 0, "Hull should be in counter-clockwise order"
  end

  def test_duplicate_points
    # Multiple points at the same location
    coords = [0.0, 0.0, 0.0, 0.0, 5.0, 0.0, 5.0, 5.0, 0.0, 5.0]
    hull = Triangles.convex_hull(coords)

    # Should handle duplicates gracefully
    assert hull.length >= 3
    assert hull.length <= 5
  end

  def test_large_point_set
    # Test with a larger set of points
    srand(123)
    points = Array.new(100) { [rand * 1000, rand * 1000] }
    coords = points.flatten

    hull = Triangles.convex_hull(coords)

    # Basic sanity checks
    assert hull.length >= 3
    assert hull.length <= 100

    # All indices should be valid
    hull.each do |idx|
      assert idx >= 0 && idx < 100
    end

    # Check for duplicates in hull
    assert_equal hull.length, hull.uniq.length, "Hull should not contain duplicate indices"
  end

  def test_tolerance_parameter
    # Test with explicit tolerance
    coords = [0.0, 0.0, 10.0, 0.0, 10.0, 10.0, 0.0, 10.0]
    hull_no_tol = Triangles.convex_hull(coords, 0.0)
    hull_with_tol = Triangles.convex_hull(coords, 0.1)

    # Both should produce valid hulls
    assert hull_no_tol.length >= 3
    assert hull_with_tol.length >= 3
  end

  def test_convex_hull_vs_delaunay_hull
    # Compare QuickHull result with Delaunay hull for validation
    srand(456)
    points = Array.new(30) { [rand * 50, rand * 50] }
    coords = points.flatten

    quickhull_result = Triangles.convex_hull(coords)
    delaunay_result = Triangles.delaunay(coords)
    delaunay_hull = delaunay_result[:hull]

    # Both should produce hulls of the same size for non-degenerate input
    # (This assumes delaunay returns correct hull)
    assert_equal quickhull_result.length, delaunay_hull.length,
                 "QuickHull and Delaunay should produce hulls of same size"
  end

  def test_negative_coordinates
    coords = [-10.0, -10.0, 10.0, -10.0, 10.0, 10.0, -10.0, 10.0, 0.0, 0.0]
    hull = Triangles.convex_hull(coords)

    # Should include all 4 corners
    assert_equal 4, hull.length
    assert_equal [0, 1, 2, 3].sort, hull.sort
  end

  def test_performance
    # Basic performance check
    srand(789)
    points = Array.new(1000) { [rand * 10000, rand * 10000] }
    coords = points.flatten

    start_time = Time.now
    hull = Triangles.convex_hull(coords)
    elapsed = Time.now - start_time

    # Should complete in reasonable time (< 1 second for 1000 points)
    assert elapsed < 1.0, "QuickHull should be fast (took #{elapsed}s)"
    assert hull.length >= 3
    puts "QuickHull processed 1000 points in #{(elapsed * 1000).round(2)}ms"
  end
end

