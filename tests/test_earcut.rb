require 'minitest/autorun'
require 'timeout'
require 'json'

require_relative '../triangles'

class EarcutTest < Minitest::Test
  FIXTURES_PATH = File.join(__dir__, 'data', 'earcut')
  EXPECTED = JSON.parse(File.read(File.join(FIXTURES_PATH, 'expected.json')))

  def test_indices_2d
    indices = Triangles.earcut([10, 0, 0, 50, 60, 60, 70, 10])
    assert_equal [1, 0, 3, 3, 2, 1], indices
  end

  def test_indices_3d
    indices = Triangles.earcut([10, 0, 0, 0, 50, 0, 60, 60, 0, 70, 10, 0], nil, 3)
    assert_equal [1, 0, 3, 3, 2, 1], indices
  end

  def test_empty
    assert_equal [], Triangles.earcut([])
  end

  # Test all fixtures from expected.json with all rotations (0, 90, 180, 270)
  EXPECTED['triangles'].each do |id, expected_triangles|
    [0, 90, 180, 270].each do |rotation|
      define_method("test_#{id.gsub('-', '_')}_rotation_#{rotation}") do
        run_fixture_test(id, expected_triangles, rotation)
      end
    end
  end

  def test_infinite_loop
    # Should not hang and should return structurally valid indices.
    indices = nil
    Timeout.timeout(5) do
      indices = Triangles.earcut([1, 2, 2, 2, 1, 2, 1, 1, 1, 2, 4, 1, 5, 1, 3, 2, 4, 2, 4, 1], [5], 2)
    end

    assert_equal 0, indices.length % 3

    n = 20 / 2
    indices.each do |idx|
      assert idx.is_a?(Integer)
      assert idx >= 0 && idx < n, "Index out of range"
    end
  end

  def test_simple_polygon_with_hole
    # Square with square hole inside
    outer = [[0, 0], [10, 0], [10, 10], [0, 10]]
    hole = [[2, 2], [2, 8], [8, 8], [8, 2]]
    
    vertices = (outer + hole).flatten
    holes = [outer.length]  # Hole starts at index 4
    
    indices = Triangles.earcut(vertices, holes, 2)
    
    # Should produce valid triangulation
    assert indices.length > 0
    assert_equal 0, indices.length % 3, "Should have complete triangles"
    
    # Verify all indices are valid
    n = (outer.length + hole.length)
    indices.each do |idx|
      assert idx >= 0 && idx < n, "Index out of range"
    end
  end

  def test_dimensions_4d
    # Test with 4D coordinates (only first 2 used)
    indices = Triangles.earcut([10, 0, 0, 0, 0, 50, 0, 0, 60, 60, 0, 0, 70, 10, 0, 0], nil, 4)
    assert_equal [1, 0, 3, 3, 2, 1], indices
  end

  def test_single_triangle_3d
    # Single triangle in 3D (z ignored)
    indices = Triangles.earcut([0, 0, 5, 10, 0, 3, 5, 10, 7], nil, 3)
    assert_equal 3, indices.length
    assert_equal [0, 1, 2].sort, indices.sort
  end

  def test_large_polygon
    # Generate a larger polygon (circle approximation)
    n = 50
    vertices = []
    n.times do |i|
      angle = 2 * Math::PI * i / n
      vertices << Math.cos(angle) * 100
      vertices << Math.sin(angle) * 100
    end
    
    indices = Triangles.earcut(vertices, nil, 2)
    
    # Should produce n-2 triangles
    expected_triangles = n - 2
    assert_equal expected_triangles * 3, indices.length
    
    # Verify deviation is small
    err = deviation(vertices, nil, 2, indices)
    assert err < 1e-10, "Deviation too large: #{err}"
  end

  def test_very_small_polygon
    # Very small polygon (scaled down)
    scale = 1e-8
    vertices = [0, 0, 1*scale, 0, 0.5*scale, 1*scale]
    indices = Triangles.earcut(vertices, nil, 2)
    
    assert_equal 3, indices.length
  end

  def test_negative_coordinates
    # Polygon with negative coordinates
    vertices = [-5, -5, -2, -5, -3.5, -2]
    indices = Triangles.earcut(vertices, nil, 2)
    
    assert_equal 3, indices.length
  end

  def test_mixed_coordinates
    # Mix of positive and negative
    vertices = [-10, -10, 10, -10, 10, 10, -10, 10]
    indices = Triangles.earcut(vertices, nil, 2)
    
    assert_equal 6, indices.length
  end

  def test_integer_pixel_coordinates
    # Common use case: integer pixel coordinates
    vertices = [0, 0, 32, 0, 32, 32, 0, 32]
    indices = Triangles.earcut(vertices, nil, 2)
    
    # Should produce 2 triangles (6 indices)
    assert_equal 6, indices.length
    
    # Both triangulations are valid:
    # [1, 0, 3, 3, 2, 1] (diagonal 0-2) or [2, 3, 0, 0, 1, 2] (diagonal 1-3)
    valid_triangulations = [
      [1, 0, 3, 3, 2, 1],  # One diagonal
      [2, 3, 0, 0, 1, 2]   # Other diagonal
    ]
    
    assert_includes valid_triangulations, indices,
      "Triangulation #{indices.inspect} is not one of the expected valid triangulations"
    
    # Verify all indices reference valid vertices
    indices.each do |idx|
      assert idx >= 0 && idx < 4, "Index #{idx} out of range"
    end
    
    # Verify the triangulation has near-zero deviation
    err = deviation(vertices, nil, 2, indices)
    assert err < 1e-10, "Deviation too large: #{err}"
  end

  private

  def run_fixture_test(id, expected_triangles, rotation = 0)
    fixture_path = File.join(FIXTURES_PATH, "#{id}.json")
    return skip("Fixture not found: #{id}") unless File.exist?(fixture_path)

    coords = JSON.parse(File.read(fixture_path))

    # Apply rotation if needed
    if rotation != 0
      theta = rotation * Math::PI / 180
      xx = Math.cos(theta).round
      xy = (-Math.sin(theta)).round
      yx = Math.sin(theta).round
      yy = Math.cos(theta).round
      coords.each do |ring|
        ring.each do |coord|
          x, y = coord[0], coord[1]
          coord[0] = xx * x + xy * y
          coord[1] = yx * x + yy * y
        end
      end
    end

    data = flatten(coords)
    indices = Triangles.earcut(data[:vertices], data[:holes], data[:dimensions])

    # Basic structural validity: complete triangles and no invalid indices.
    assert_equal 0, indices.length % 3
    n = data[:vertices].length / data[:dimensions]
    indices.each do |idx|
      assert idx.is_a?(Integer)
      assert idx >= 0 && idx < n, "Index #{idx} out of range"
    end

    err = deviation(data[:vertices], data[:holes], data[:dimensions], indices)

    # Get expected deviation - use rotation-specific if available, otherwise default
    expected_deviation = if rotation != 0 && EXPECTED['errors-with-rotation'] && EXPECTED['errors-with-rotation'][id]
                           EXPECTED['errors-with-rotation'][id]
                         else
                           EXPECTED['errors'][id] || 0
                         end

    num_triangles = indices.length / 3

    # Only check triangle count for rotation 0
    if rotation == 0
      assert_equal expected_triangles, num_triangles, "#{id}: #{num_triangles} triangles when expected #{expected_triangles}"
    end

    # Check deviation for all rotations if we expect triangles
    if expected_triangles > 0
      assert err <= expected_deviation, "#{id} rotation #{rotation}: deviation #{err} > #{expected_deviation}"
    else
      # When we expect no triangles, ensure it truly returns nothing.
      assert_equal [], indices
    end
  end

  # Flatten nested polygon coordinates into a flat array
  def flatten(data)
    vertices = []
    holes = []
    dimensions = data[0][0].length
    hole_index = 0
    prev_len = 0

    data.each do |ring|
      ring.each do |p|
        d = 0
        while d < dimensions
          vertices << p[d]
          d += 1
        end
      end
      if prev_len > 0
        hole_index += prev_len
        holes << hole_index
      end
      prev_len = ring.length
    end

    { vertices: vertices, holes: holes, dimensions: dimensions }
  end

  # Calculate deviation between polygon area and triangulation area
  def deviation(data, hole_indices, dim, triangles)
    has_holes = hole_indices && hole_indices.length > 0
    outer_len = has_holes ? hole_indices[0] * dim : data.length

    polygon_area = signed_area(data, 0, outer_len, dim).abs
    if has_holes
      i = 0
      len = hole_indices.length
      while i < len
        start = hole_indices[i] * dim
        end_idx = i < len - 1 ? hole_indices[i + 1] * dim : data.length
        polygon_area -= signed_area(data, start, end_idx, dim).abs
        i += 1
      end
    end

    triangles_area = 0.0
    i = 0
    while i < triangles.length
      a = triangles[i] * dim
      b = triangles[i + 1] * dim
      c = triangles[i + 2] * dim
      triangles_area += ((data[a] - data[c]) * (data[b + 1] - data[a + 1]) -
                         (data[a] - data[b]) * (data[c + 1] - data[a + 1])).abs
      i += 3
    end

    return 0 if polygon_area == 0 && triangles_area == 0
    ((triangles_area - polygon_area) / polygon_area).abs
  end

  def signed_area(data, start, end_idx, dim)
    sum = 0.0
    j = end_idx - dim
    i = start
    while i < end_idx
      sum += (data[j] - data[i]) * (data[i + 1] + data[j + 1])
      j = i
      i += dim
    end
    sum
  end
end
