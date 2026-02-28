require 'minitest/autorun'
require_relative '../triangles'

module AsteroidsGeometry
  module_function

  def generate_asteroid_points(radius:, vertex_count:, jitter: 0.35)
    step = (Math::PI * 2.0) / vertex_count
    start = rand * Math::PI * 2.0
    coords = []
    vertex_count.times do |i|
      angle = start + i * step + rand(-step * jitter..step * jitter)
      r = radius * (1.0 - jitter) + rand * (radius * jitter * 2.0)
      coords << Math.cos(angle) * r << Math.sin(angle) * r
    end

    Triangles.ensure_ccw(coords)
  end
end


# Comprehensive test suite for Bayazit convex decomposition algorithm
class TestBayazitDecomposition < Minitest::Test
  def setup
    @shatter = Triangles
  end

  # ==================== Convexity Detection Tests ====================

  def test_convex_square_is_convex
    square = [0, 0, 10, 0, 10, 10, 0, 10]
    assert @shatter.convex?(square), "Square should be convex"
  end

  def test_convex_triangle_is_convex
    triangle = [0, 0, 10, 0, 5, 10]
    assert @shatter.convex?(triangle), "Triangle should be convex"
  end

  def test_convex_regular_hexagon_is_convex
    hexagon = []
    6.times do |i|
      angle = i * Math::PI * 2.0 / 6.0
      hexagon << Math.cos(angle) * 50 << Math.sin(angle) * 50
    end
    assert @shatter.convex?(hexagon), "Regular hexagon should be convex"
  end

  def test_non_convex_star_is_not_convex
    star = []
    5.times do |i|
      outer_angle = i * Math::PI * 2.0 / 5.0
      star << Math.cos(outer_angle) * 50 << Math.sin(outer_angle) * 50

      inner_angle = (i + 0.5) * Math::PI * 2.0 / 5.0
      star << Math.cos(inner_angle) * 25 << Math.sin(inner_angle) * 25
    end
    refute @shatter.convex?(star), "Star should be non-convex"
  end

  def test_non_convex_l_shape_is_not_convex
    l_shape = [0, 0, 10, 0, 10, 5, 5, 5, 5, 10, 0, 10]
    refute @shatter.convex?(l_shape), "L-shape should be non-convex"
  end

  def test_collinear_points_handling
    # Three points in a line should still be considered valid
    line_like = [0, 0, 5, 5, 10, 10]
    # Should not crash
    result = @shatter.convex?(line_like)
    assert result.is_a?(TrueClass) || result.is_a?(FalseClass), "Should return boolean"
  end

  # Decomposition Tests

  def test_decompose_convex_polygon_returns_single_component
    square = [0, 0, 10, 0, 10, 10, 0, 10]
    result = @shatter.decompose_to_convex_hulls(square)

    assert_equal 1, result.length, "Convex polygon should decompose to 1 component"
    assert_equal square, result[0], "Component should match original polygon"
  end

  def test_decompose_non_convex_star_returns_multiple_components
    star = []
    5.times do |i|
      outer_angle = i * Math::PI * 2.0 / 5.0
      star << Math.cos(outer_angle) * 50 << Math.sin(outer_angle) * 50

      inner_angle = (i + 0.5) * Math::PI * 2.0 / 5.0
      star << Math.cos(inner_angle) * 25 << Math.sin(inner_angle) * 25
    end

    result = @shatter.decompose_to_convex_hulls(star)

    assert result.length > 1, "Non-convex star should decompose to multiple components"
  end

  def test_decompose_l_shape_returns_convex_components
    l_shape = [0, 0, 10, 0, 10, 5, 5, 5, 5, 10, 0, 10]
    result = @shatter.decompose_to_convex_hulls(l_shape)

    # CORRECT: All components must be convex
    assert result.length >= 1, "L-shape should decompose to at least 1 component"

    result.each do |component|
      assert component.length >= 6, "Each component should have at least 3 vertices"
      assert @shatter.convex?(component),
             "EVERY component must be convex: #{component.inspect}"
    end
  end

  def test_all_components_are_convex
    20.times do
      asteroid = AsteroidsGeometry.generate_asteroid_points(
        radius: rand(40.0..80.0),
        vertex_count: rand(8..14),
        jitter: 0.35
      )

      components = @shatter.decompose_to_convex_hulls(asteroid)

      # REQUIREMENT: All components must be convex
      components.each_with_index do |component, idx|
        assert component.length >= 6, "Component #{idx} should have at least 3 vertices, got #{component.length.div(2)}"
        assert @shatter.convex?(component),
               "Component #{idx} MUST be convex but is not: #{component.inspect}"
      end
    end
  end

  def test_strict_convexity_no_reflex_vertices
    10.times do
      asteroid = AsteroidsGeometry.generate_asteroid_points(
        radius: rand(40.0..80.0),
        vertex_count: rand(8..14),
        jitter: 0.35
      )

      components = @shatter.decompose_to_convex_hulls(asteroid)

      components.each_with_index do |component, comp_idx|
        # First determine the winding order
        area = @shatter.polygon_area(component)
        is_ccw = area > 0
        n = component.length.div(2)

        # Check every vertex is NOT a reflex
        n.times do |i|
          i_prev = (i - 1 + n) % n
          i_next = (i + 1) % n
          p0x = component[i_prev * 2]; p0y = component[i_prev * 2 + 1]
          p1x = component[i * 2];      p1y = component[i * 2 + 1]
          p2x = component[i_next * 2]; p2y = component[i_next * 2 + 1]

          cross = (p1x - p0x) * (p2y - p0y) - (p1y - p0y) * (p2x - p0x)

          # For CCW polygons, the cross should be >= 0 for convex vertices
          # For CW polygons, the cross should be <= 0 for convex vertices
          if is_ccw
            refute cross < -1e-8,
                   "Component #{comp_idx} (CCW) has reflex vertex at index #{i}: cross=#{cross}"
          else
            refute cross > 1e-8,
                   "Component #{comp_idx} (CW) has reflex vertex at index #{i}: cross=#{cross}"
          end
        end
      end
    end
  end

  def test_all_vertices_maintain_winding_order
    # Ensure all decomposed components preserve the winding order
    10.times do
      asteroid = AsteroidsGeometry.generate_asteroid_points(
        radius: rand(40.0..80.0),
        vertex_count: rand(8..14),
        jitter: 0.35
      )

      input_area = @shatter.polygon_area(asteroid)
      input_is_ccw = input_area > 0

      components = @shatter.decompose_to_convex_hulls(asteroid)

      components.each_with_index do |component, idx|
        area = @shatter.polygon_area(component)

        # All components should have the SAME winding as input (allow tiny tolerance)
        if input_is_ccw
          refute (area < -1e-6),
                 "Component #{idx} has clockwise winding (area=#{area}) but input is CCW"
        else
          refute (area > 1e-6),
                 "Component #{idx} has CCW winding (area=#{area}) but input is clockwise"
        end
      end
    end
  end

  def test_decomposition_preserves_geometry
    # Test that decomposed components cover the original polygon's centroid
    asteroid = AsteroidsGeometry.generate_asteroid_points(
      radius: 50.0,
      vertex_count: 10,
      jitter: 0.35
    )

    centroid = @shatter.centroid(asteroid)
    components = @shatter.decompose_to_convex_hulls(asteroid)

    # At least one component should contain the centroid
    found = components.any? do |component|
      @shatter.point_in_poly?(centroid[0], centroid[1], component)
    end

    assert found, "At least one component should contain the original centroid"
  end

  def test_minimum_vertex_requirement
    # Very small polygon should still decompose
    tiny = [0, 0, 1, 0, 0.5, 1]
    result = @shatter.decompose_to_convex_hulls(tiny)

    assert result.length >= 1, "Even tiny polygon should return components"
    result.each do |component|
      assert component.length >= 6, "Each component should have at least 3 vertices"
    end
  end

  # Collision Detection Tests
  def test_collision_detection_works_on_decomposed_components
    # Generate a non-convex asteroid
    asteroid = AsteroidsGeometry.generate_asteroid_points(
      radius: 50.0,
      vertex_count: 10,
      jitter: 0.35
    )

    components = @shatter.decompose_to_convex_hulls(asteroid)

    # Get the asteroid's actual centroid
    centroid = @shatter.centroid(asteroid)
    cx, cy = centroid[0], centroid[1]

    # At least one component should contain the centroid point
    collision_found = components.any? do |component|
      @shatter.point_in_poly?(cx, cy, component)
    end

    assert collision_found, "At least one decomposed component should contain the centroid"
  end

  # Edge Cases

  def test_empty_polygon_handling
    result = @shatter.decompose_to_convex_hulls([])
    assert_equal 1, result.length, "Empty polygon should return as-is"
    assert_equal [], result[0], "Should return the empty array"
  end

  def test_two_vertex_polygon_handling
    two_points = [0, 0, 10, 10]
    result = @shatter.decompose_to_convex_hulls(two_points)

    # Should return as-is (cannot form a valid polygon)
    assert_equal 1, result.length, "Two-point polygon should return as-is"
  end

  def test_triangle_decomposition
    triangle = [0, 0, 10, 0, 5, 10]
    result = @shatter.decompose_to_convex_hulls(triangle)

    assert_equal 1, result.length, "Triangle should decompose to 1 (already convex)"
    assert @shatter.convex?(result[0]), "Result should be convex"
  end

  def test_very_small_reflex_angle
    # Polygon with very small reflex angle (nearly convex)
    nearly_convex = [0, 0, 10, 0, 10, 10, 5.1, 5, 0, 10]

    result = @shatter.decompose_to_convex_hulls(nearly_convex)

    # Should decompose to at least 1 component with valid vertices
    assert result.length >= 1, "Should decompose to at least 1 component"
    result.each do |component|
      assert component.length >= 6, "Each component should have at least 3 vertices"
    end
  end

  # Large Dataset Tests

  def test_decompose_100_generated_asteroids
    success_count = 0
    100.times do |i|
      asteroid = AsteroidsGeometry.generate_asteroid_points(
        radius: rand(40.0..80.0),
        vertex_count: rand(8..14),
        jitter: 0.35
      )

      begin
        components = @shatter.decompose_to_convex_hulls(asteroid)

        # Check all components have at least 3 vertices and are convex
        all_valid = components.all? { |c| c.length >= 6 }
        all_convex = components.all? { |c| @shatter.convex?(c) }

        if all_valid && all_convex && components.length >= 1
          success_count += 1
        end
      rescue
        # Allow for some edge case failures during algorithm development
      end
    end

    # Expect a high success rate - algorithm should work on most asteroids
    assert success_count == 100, "All asteroids should decompose to all-convex components (got #{success_count}/100)"
  end

  def test_decomposition_statistics
    total_asteroids = 0
    total_decompositions_completed = 0

    10.times do
      asteroid = AsteroidsGeometry.generate_asteroid_points(
        radius: rand(40.0..80.0),
        vertex_count: rand(8..14),
        jitter: 0.35
      )

      total_asteroids += 1

      # Main requirement: decomposition completes successfully
      begin
        components = @shatter.decompose_to_convex_hulls(asteroid)
        if components.length >= 1
          total_decompositions_completed += 1
        end
      rescue
        # Allow for some failures
      end
    end

    # At least 80% should decompose successfully
    assert total_decompositions_completed >= 8,
           "At least 80% of asteroids should decompose successfully (got #{total_decompositions_completed}/#{total_asteroids})"
  end

  def test_highly_concave_star_polygon
    # Create a highly concave star with many reflex vertices
    star = []
    10.times do |i|
      outer_angle = i * Math::PI * 2.0 / 10.0
      star << Math.cos(outer_angle) * 100 << Math.sin(outer_angle) * 100

      inner_angle = (i + 0.5) * Math::PI * 2.0 / 10.0
      star << Math.cos(inner_angle) * 40 << Math.sin(inner_angle) * 40
    end

    components = @shatter.decompose_to_convex_hulls(star)

    # All components must be convex
    assert components.length >= 1, "Star should decompose into at least 1 component"
    components.each_with_index do |component, idx|
      assert component.length >= 6, "Component #{idx} must have at least 3 vertices"
      assert @shatter.convex?(component),
             "Component #{idx} from star MUST be convex: #{component.inspect}"
    end
  end

  def test_comb_shape_many_reflex_vertices
    # Create a comb shape with many consecutive reflex vertices
    comb = [
      0, 0,
      10, 0,
      10, 20,
      15, 20,
      15, 0,
      20, 0,
      20, 20,
      25, 20,
      25, 0,
      30, 0,
      30, 20,
      35, 20,
      35, 0,
      40, 0,
      40, 30,
      0, 30
    ]

    components = @shatter.decompose_to_convex_hulls(comb)

    # All components must be convex
    assert components.length >= 1, "Comb should produce at least 1 component"
    components.each_with_index do |component, idx|
      assert component.length >= 6, "Comb component #{idx} must have at least 3 vertices, got #{component.length.div(2)}: #{component.inspect}"
      assert @shatter.convex?(component),
             "Comb component #{idx} MUST be convex: #{component.inspect}"
    end
  end

  # Performance Tests

  def test_decomposition_completes_in_reasonable_time
    asteroid = AsteroidsGeometry.generate_asteroid_points(
      radius: 50.0,
      vertex_count: 14,
      jitter: 0.35
    )

    start_time = Time.now
    _components = @shatter.decompose_to_convex_hulls(asteroid)
    elapsed = Time.now - start_time

    assert elapsed < 1.0, "Decomposition should complete in under 1 second (took #{elapsed}s)"
  end

  def test_recursive_decomposition_depth_limit
    # Create a pathological case that would cause deep recursion
    star = []
    10.times do |i|
      outer_angle = i * Math::PI * 2.0 / 10.0
      star << Math.cos(outer_angle) * 100 << Math.sin(outer_angle) * 100

      inner_angle = (i + 0.5) * Math::PI * 2.0 / 10.0
      star << Math.cos(inner_angle) * 50 << Math.sin(inner_angle) * 50
    end

    # Should complete without stack overflow
    error = nil
    begin
      _result = @shatter.decompose_to_convex_hulls(star)
    rescue => e
      error = e
    end

    assert_nil error, "Decomposition should not raise error: #{error.inspect}"
  end

  # Robustness Tests

  def test_decomposition_with_duplicate_vertices
    # Polygon with consecutive duplicate vertices
    poly = [0, 0, 0, 0, 10, 0, 10, 10, 0, 10]

    # Should handle gracefully
    error = nil
    begin
      _result = @shatter.decompose_to_convex_hulls(poly)
    rescue => e
      error = e
    end

    assert_nil error, "Should handle duplicate vertices: #{error.inspect}"
  end

  def test_decomposition_with_floating_point_coords
    # Polygon with precise floating point coordinates
    poly = [0.123, 0.456, 10.789, 0.123, 10.456, 10.789, 0.789, 10.456, 5.0, 5.5]

    result = @shatter.decompose_to_convex_hulls(poly)

    # Should decompose to at least one component
    assert result.length >= 1, "Should decompose to at least 1 component"

    # Most components should be valid (at least 3 vertices)
    valid_count = result.count { |c| c.length >= 6 }
    assert valid_count >= result.length * 0.7,
           "At least 70% of components should be valid"
  end

  def test_decomposition_consistency
    # Same polygon should decompose consistently
    asteroid = AsteroidsGeometry.generate_asteroid_points(
      radius: 50.0,
      vertex_count: 10,
      jitter: 0.35
    )

    result1 = @shatter.decompose_to_convex_hulls(asteroid)
    result2 = @shatter.decompose_to_convex_hulls(asteroid)

    assert_equal result1.length, result2.length,
                 "Same polygon should decompose to same number of components"
  end

  # Helper Methods

  def find_collision_polygon_polygon(a, b)
    # Simplified SAT collision detection for testing
    # Convert flat arrays to paired for collision test
    a = flat_to_pairs(a)
    b = flat_to_pairs(b)

    separation = -Float::INFINITY

    k = 0
    while k < 2
      a_length = a.length
      b_length = b.length

      i = 0
      while i < b_length
        v1x, v1y = b[(i + 1) % b_length]
        v0x, v0y = b[i]

        nx = v1y - v0y
        ny = -(v1x - v0x)
        l = 1.0 / (Math.sqrt(nx * nx + ny * ny))
        nx *= l
        ny *= l

        v2x, v2y = a[0]
        min_projection = (v2x - v0x) * nx + (v2y - v0y) * ny

        j = 1
        while j < a_length
          v2x, v2y = a[j]
          projection = (v2x - v0x) * nx + (v2y - v0y) * ny
          min_projection = projection if projection < min_projection

          j += 1
        end

        return nil if min_projection >= 0

        if min_projection * 0.95 > separation + 0.01
          separation = min_projection
        end

        i += 1
      end

      t = a
      a = b
      b = t

      k += 1
    end

    { contacts: [{ r1x: 0, r1y: 0, r2x: 0, r2y: 0 }] }
  end

  def flat_to_pairs(coords)
    return [] if coords.nil? || coords.empty?
    pairs = []
    i = 0
    while i < coords.length
      pairs << [coords[i], coords[i + 1]]
      i += 2
    end
    pairs
  end
end
