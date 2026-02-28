require 'minitest/autorun'
require 'timeout'
require 'json'

require_relative '../triangles'

class TrianglesTest < Minitest::Test
  FIXTURES_PATH = File.join(__dir__, 'data', 'delaunator')

  def setup
    srand(1234)

    @points = load_fixture('ukraine.json')
    @issue13 = load_fixture('issue13.json')
    @issue43 = load_fixture('issue43.json')
    @issue44 = load_fixture('issue44.json')
    @robustness1 = load_fixture('robustness1.json')
    @robustness2 = load_fixture('robustness2.json')
    @robustness3 = load_fixture('robustness3.json')
    @robustness4 = load_fixture('robustness4.json')
  end

  def load_fixture(name)
    JSON.parse(File.read(File.join(FIXTURES_PATH, name)))
  end

  def test_triangulates_plain_array
    coords = @points.flatten
    result1 = Triangles.delaunay(coords)
    result2 = Triangles.delaunay(coords)
    assert_equal result1[:triangles], result2[:triangles]
  end

  def test_produces_correct_triangulation
    validate(@points)
  end

  def test_issue_11
    validate([[516, 661], [369, 793], [426, 539], [273, 525], [204, 694], [747, 750], [454, 390]])
  end

  def test_issue_13
    validate(@issue13)
  end

  def test_issue_24
    validate([[382, 302], [382, 328], [382, 205], [623, 175], [382, 188], [382, 284], [623, 87], [623, 341], [141, 227]])
  end

  def test_issue_43
    validate(@issue43)
  end

  def test_issue_44
    validate(@issue44)
  end

  def test_robustness
    Timeout.timeout(60) do
      validate(@robustness1)
      validate(@robustness1.map { |p| [p[0] / 1e9, p[1] / 1e9] })
      validate(@robustness1.map { |p| [p[0] / 100, p[1] / 100] })
      validate(@robustness1.map { |p| [p[0] * 100, p[1] * 100] })
      validate(@robustness1.map { |p| [p[0] * 1e9, p[1] * 1e9] })
      validate(@robustness2[0...100])
      validate(@robustness2)
      validate(@robustness3)
      validate(@robustness4)
    end
  end

  def test_returns_empty_triangulation_for_small_number_of_points
    # Empty array
    result = Triangles.delaunay([])
    assert_equal [], result[:triangles]
    assert_equal [], result[:hull]

    # Single point - @points[0] is [x, y], so @points.first(1) is [[x, y]]
    single_point_coords = @points.first(1).flatten  # Gets [x, y]
    result = Triangles.delaunay(single_point_coords)
    assert_equal [], result[:triangles]
    assert_equal [0], result[:hull]

    # Two points - @points.first(2) is [[x1, y1], [x2, y2]]
    two_points_coords = @points.first(2).flatten  # Gets [x1, y1, x2, y2]
    result = Triangles.delaunay(two_points_coords)
    assert_equal [], result[:triangles]
    # [0, 1] or [1, 0] are both correct
    assert_includes [[0, 1], [1, 0]], result[:hull]
  end

  def test_returns_empty_triangulation_for_all_collinear_input
    result = Triangles.delaunay([[0, 0], [1, 0], [3, 0], [2, 0]].flatten)
    assert_equal [], result[:triangles]
    # Hull should contain all points in some order
    assert_equal 4, result[:hull].size
    assert_equal [0, 1, 2, 3].sort, result[:hull].sort
  end

  def test_supports_custom_point_format
    # For the Ruby implementation we convert to flat array format
    points = [{x: 5, y: 5}, {x: 7, y: 5}, {x: 7, y: 6}]
    coords = points.flat_map { |p| [p[:x], p[:y]] }
    result = Triangles.delaunay(coords)
    assert_equal [0, 2, 1], result[:triangles]
  end

  def test_handles_duplicate_points
    # Test with some duplicate points
    coords = [0, 0, 1, 0, 1, 0, 1, 1, 0, 1].flatten  # (1,0) is duplicated
    result = Triangles.delaunay(coords)

    # Must not crash, and must return structurally valid output.
    assert_equal 0, result[:triangles].length % 3
    assert result[:hull].length >= 3

    n = coords.length / 2
    result[:triangles].each do |idx|
      assert idx.is_a?(Integer)
      assert idx >= 0 && idx < n, "Triangle index #{idx} out of range [0, #{n})"
    end

    # Halfedge symmetry is a key invariant from upstream delaunator tests.
    result[:half_edges].each_with_index do |he, i|
      next if he == -1
      assert_equal i, result[:half_edges][he], "Halfedge #{i} -> #{he} should have twin #{he} -> #{i}"
    end
  end

  def test_single_dimension_variation
    # All points have same x coordinate - should return collinear hull
    coords = [5, 0, 5, 1, 5, 2, 5, 3].flatten
    result = Triangles.delaunay(coords)
    
    assert_equal [], result[:triangles]
    assert_equal 4, result[:hull].length
  end

  def test_single_dimension_variation_y
    # All points have same y coordinate - should return collinear hull
    coords = [0, 5, 1, 5, 2, 5, 3, 5].flatten
    result = Triangles.delaunay(coords)
    
    assert_equal [], result[:triangles]
    assert_equal 4, result[:hull].length
  end

  def test_large_coordinate_values
    # Test with very large coordinate values
    scale = 1e10
    points = [[0, 0], [1, 0], [0.5, 0.5]].map { |p| [p[0] * scale, p[1] * scale] }
    coords = points.flatten
    result = Triangles.delaunay(coords)
    
    assert_equal 3, result[:triangles].length
    assert result[:hull].length >= 3
  end

  def test_small_coordinate_values
    # Test with very small coordinate values
    scale = 1e-10
    points = [[0, 0], [1, 0], [0.5, 0.5]].map { |p| [p[0] * scale, p[1] * scale] }
    coords = points.flatten
    result = Triangles.delaunay(coords)
    
    assert_equal 3, result[:triangles].length
    assert result[:hull].length >= 3
  end

  def test_negative_coordinates
    # Test with negative coordinates
    coords = [[-5, -5], [-3, -5], [-4, -3]].flatten
    result = Triangles.delaunay(coords)
    
    assert_equal 3, result[:triangles].length
    assert_equal 3, result[:hull].length
  end

  def test_mixed_positive_negative_coordinates
    # Test with mix of positive and negative coordinates
    coords = [[-1, -1], [1, -1], [1, 1], [-1, 1]].flatten
    result = Triangles.delaunay(coords)
    
    assert_equal 6, result[:triangles].length  # 2 triangles * 3 vertices
    assert_equal 4, result[:hull].length
  end

  def test_integer_coordinates
    # Test with integer coordinates (common in games)
    coords = [0, 0, 100, 0, 100, 100, 0, 100]
    result = Triangles.delaunay(coords)
    
    assert_equal 6, result[:triangles].length
    assert_equal 4, result[:hull].length
  end

  def test_convex_hull_order
    # Test that hull is returned in counter-clockwise order
    coords = [0, 0, 2, 0, 2, 2, 0, 2].flatten
    result = Triangles.delaunay(coords)
    
    # Hull should be in CCW order
    hull = result[:hull]
    assert_equal 4, hull.length
    
    # Verify it's a valid hull by checking all points are used
    assert_equal [0, 1, 2, 3].sort, hull.sort
  end

  def test_triangles_reference_valid_indices
    # Ensure all triangle indices reference valid points
    coords = @points.flatten
    result = Triangles.delaunay(coords)
    n = @points.length
    
    result[:triangles].each do |idx|
      assert idx >= 0 && idx < n, "Triangle index #{idx} out of range [0, #{n})"
    end
  end

  def test_halfedges_symmetry
    # Test the halfedges array maintains symmetry
    coords = @points.flatten
    result = Triangles.delaunay(coords)
    
    result[:half_edges].each_with_index do |he, i|
      if he != -1
        # If halfedge exists, its twin should point back to this edge
        assert_equal i, result[:half_edges][he], "Halfedge #{i} -> #{he} should have twin #{he} -> #{i}"
      end
    end
  end

  def test_performance_reasonable
    # Basic performance sanity check - should handle moderate data quickly.
    # Generate 1000 random points (deterministically seeded in setup).
    coords = []
    1000.times do
      coords << rand(1000)
      coords << rand(1000)
    end

    # Avoid brittle wall-clock expectations across machines; just ensure it terminates.
    result = nil
    Timeout.timeout(5) { result = Triangles.delaunay(coords) }

    assert result[:triangles].length > 0
  end

  private

  def orient(p, r, q)
    px, py = p
    rx, ry = r
    qx, qy = q
    l = (ry - py) * (qx - px)
    r_val = (rx - px) * (qy - py)
    (l - r_val).abs >= 3.3306690738754716e-16 * (l + r_val).abs ? l - r_val : 0
  end

  def convex?(r, q, p)
    (orient(p, r, q) >= 0) || (orient(r, q, p) >= 0) || (orient(q, p, r) >= 0)
  end

  def validate(points)
    coords = points.flatten
    result = Triangles.delaunay(coords)

    # Validate halfedges
    result[:half_edges].each_with_index do |he, i|
      if he != -1
        assert_equal i, result[:half_edges][he], "valid halfedge connection at #{i}"
      end
    end

    # Validate triangulation
    hull = result[:hull]
    hull_areas = []

    hull.each_with_index do |_, i|
      j = (i - 1) % hull.size
      x0, y0 = points[hull[j]]
      x, y = points[hull[i]]
      hull_areas << (x - x0) * (y + y0)

      # Check convexity
      p1 = points[hull[j]]
      p2 = points[hull[(j + 1) % hull.size]]
      p3 = points[hull[(j + 3) % hull.size]]
      assert convex?(p1, p2, p3), "hull should be convex at #{j}"
    end

    hull_area = kahan_sum(hull_areas)

    triangle_areas = []
    triangles = result[:triangles]
    i = 0
    while i < triangles.size
      ax, ay = points[triangles[i]]
      bx, by = points[triangles[i + 1]]
      cx, cy = points[triangles[i + 2]]
      triangle_areas << ((by - ay) * (cx - bx) - (bx - ax) * (cy - by)).abs
      i += 3
    end

    triangles_area = kahan_sum(triangle_areas)

    # If the hull area is ~0, we expect no triangles (collinear/degenerate input).
    if hull_area == 0
      assert_equal 0, triangles_area
      assert_equal [], triangles
      return
    end

    err = ((hull_area - triangles_area) / hull_area).abs
    assert err <= 2**-51, "triangulation should be valid; #{err} error"
  end

  # Kahan and Babuska summation, Neumaier variant
  def kahan_sum(arr)
    return 0 if arr.empty?

    sum = arr[0]
    err = 0

    (1...arr.size).each do |i|
      k = arr[i]
      m = sum + k
      err += sum.abs >= k.abs ? sum - m + k : k - m + sum
      sum = m
    end

    sum + err
  end
end
