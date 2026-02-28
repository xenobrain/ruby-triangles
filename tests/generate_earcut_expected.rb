#!/usr/bin/env ruby
# Generates tests/data/earcut/expected.json from the actual earcut output.
# Run once: ruby tests/generate_earcut_expected.rb

require_relative '../triangles'
require 'json'

FIXTURES_PATH = File.join(__dir__, 'data', 'earcut')

def flatten_polygon(data)
  vertices = []
  dimensions = data[0][0].length
  holes = []
  hole_index = 0
  prev_len = 0

  data.each do |ring|
    ring.each { |p| p.each { |v| vertices << v } }
    if prev_len > 0
      hole_index += prev_len
      holes << hole_index
    end
    prev_len = ring.length
  end

  { vertices: vertices, holes: holes.empty? ? nil : holes, dimensions: dimensions }
end

def deviation(data, hole_indices, dim, triangles)
  has_holes = hole_indices && !hole_indices.empty?
  outer_len = has_holes ? hole_indices[0] * dim : data.length
  polygon_area = signed_area(data, 0, outer_len, dim).abs

  if has_holes
    hole_indices.each_with_index do |start_idx, i|
      start = start_idx * dim
      end_idx = i < hole_indices.length - 1 ? hole_indices[i + 1] * dim : data.length
      polygon_area -= signed_area(data, start, end_idx, dim).abs
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

  return 0.0 if polygon_area == 0 && triangles_area == 0
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

def rotate_coords(data, degrees)
  theta = degrees * Math::PI / 180
  xx = Math.cos(theta).round
  xy = (-Math.sin(theta)).round
  yx = Math.sin(theta).round
  yy = Math.cos(theta).round
  data.map do |ring|
    ring.map do |coord|
      x, y = coord[0], coord[1]
      [xx * x + xy * y, yx * x + yy * y]
    end
  end
end

triangles_map = {}
errors_map = {}

Dir.glob(File.join(FIXTURES_PATH, '*.json')).sort.each do |path|
  name = File.basename(path, '.json')
  next if name == 'expected'

  begin
    data = JSON.parse(File.read(path))
    flat = flatten_polygon(data)
    indices = Triangles.earcut(flat[:vertices], flat[:holes], flat[:dimensions])
    num_triangles = indices.length / 3
    triangles_map[name] = num_triangles

    # Measure max error across all rotations
    max_err = 0.0
    [0, 90, 180, 270].each do |rotation|
      rotated = rotation == 0 ? data : rotate_coords(data, rotation)
      rflat = flatten_polygon(rotated)
      ridx = Triangles.earcut(rflat[:vertices], rflat[:holes], rflat[:dimensions])
      rerr = deviation(rflat[:vertices], rflat[:holes], rflat[:dimensions], ridx)
      max_err = rerr if rerr > max_err
    end

    # Set expected error: exactly 0 if all rotations produce exact result,
    # otherwise give 2x headroom with a minimum floor.
    if max_err == 0.0
      errors_map[name] = 0
    else
      errors_map[name] = [max_err * 2.0, 1e-14].max
    end

    puts "%-40s %4d triangles  max_err=%.2e  expected=%.2e" % [name, num_triangles, max_err, errors_map[name]]
  rescue => e
    puts "ERROR #{name}: #{e.message}"
    puts e.backtrace.first(3).join("\n")
  end
end

result = {
  "triangles" => triangles_map.sort.to_h,
  "errors"    => errors_map.sort.to_h
}

out = File.join(FIXTURES_PATH, 'expected.json')
File.write(out, JSON.pretty_generate(result))
puts "\nWrote #{out}"
