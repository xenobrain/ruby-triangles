require 'minitest/autorun'
require 'json'
require 'set'

require_relative '../triangles'

class TestConstrain < Minitest::Test
  FIXTURES_PATH = File.join(__dir__, 'data', 'cdt')

  # Load all test files
  def self.load_test_files
    files = []
    Dir.glob(File.join(FIXTURES_PATH, '**/*.json')).each do |path|
      data = JSON.parse(File.read(path))
      name = File.basename(path, '.json')
      files << {
        name: name,
        path: path,
        points: data['points'],
        edges: data['edges'],
        error: data['error']
      }
    end
    files
  end

  TEST_FILES = load_test_files

  # Helper: next edge in triangle
  def next_edge(e)
    (e % 3 == 2) ? e - 2 : e + 1
  end

  # Helper: previous edge in triangle
  def prev_edge(e)
    (e % 3 == 0) ? e + 2 : e - 1
  end

  # Validate that half-edges link back correctly
  def validate_halfedges(del)
    del[:triangles].length.times do |edg|
      adj = del[:half_edges][edg]
      next if adj == -1

      assert_equal edg, del[:half_edges][adj],
        "Half-edge #{edg} -> #{adj} should link back"

      e1 = del[:triangles][edg]
      e2 = del[:triangles][next_edge(edg)]
      a1 = del[:triangles][adj]
      a2 = del[:triangles][next_edge(adj)]

      assert(e1 == a2 || e2 == a1,
        "Half-edges #{edg}/#{adj} should share endpoints: (#{e1},#{e2}) vs (#{a1},#{a2})")
    end
  end

  # Validate vertex map - every vertex should be reachable
  def validate_vert_map(con, points)
    del = con[:del]
    num_points = points.length

    # Build map of which edges point to each vertex
    edge_map = Hash.new { |h, k| h[k] = Set.new }
    del[:triangles].length.times do |edg|
      p1 = del[:triangles][edg]
      edge_map[p1] << prev_edge(edg)
    end

    num_points.times do |i|
      incoming = edge_map[i]
      next if incoming.empty?

      start = con[:vert_map][i]
      assert start, "Vertex #{i} should have an entry in vert_map"

      # Walk around the vertex
      edg = start
      visited = Set.new
      loop do
        assert incoming.include?(edg),
          "Edge #{edg} should be incoming to vertex #{i}"
        visited << edg

        nxt = next_edge(edg)
        adj = del[:half_edges][nxt]
        edg = adj
        break if edg == -1 || edg == start
      end

      # All incoming edges should be visited
      assert_equal incoming, visited,
        "All edges to vertex #{i} should be visited: expected #{incoming.to_a}, got #{visited.to_a}"
    end
  end

  # Validate that flips are cleared (no pending flips)
  def validate_flips_cleared(con)
    assert con[:flips].empty?, "Flips should be cleared after constraining"
  end

  # Validate that constrained edges are consistent between half-edges
  def validate_consd_consistency(con)
    del = con[:del]
    del[:triangles].length.times do |edg|
      adj = del[:half_edges][edg]
      next if adj == -1

      assert_equal Triangles::TriVis.cdt_edge_constrained?(con, edg),
                   Triangles::TriVis.cdt_edge_constrained?(con, adj),
        "Constrained status should be consistent for half-edges #{edg}/#{adj}"
    end
  end

  # Check if two segments intersect (robust)
  # Note: For collinear segments, we only report intersection if they actually cross,
  # not if they just overlap/contain each other on the same line
  def segments_intersect?(p1, p2, p3, p4)
    x1, y1 = p1
    x2, y2 = p2
    x3, y3 = p3
    x4, y4 = p4

    o1 = Triangles::Geometry.triangle_signed_area(x1, y1, x2, y2, x3, y3)
    o2 = Triangles::Geometry.triangle_signed_area(x1, y1, x2, y2, x4, y4)
    o3 = Triangles::Geometry.triangle_signed_area(x3, y3, x4, y4, x1, y1)
    o4 = Triangles::Geometry.triangle_signed_area(x3, y3, x4, y4, x2, y2)

    # Use tolerance for floating point comparison
    eps = 1e-10
    o1_zero = o1.abs < eps
    o2_zero = o2.abs < eps
    o3_zero = o3.abs < eps
    o4_zero = o4.abs < eps

    # If all four points are collinear, this is NOT a crossing intersection
    # (collinear overlapping segments are allowed in CDT - they just share the same line)
    if o1_zero && o2_zero && o3_zero && o4_zero
      return false  # Collinear segments don't "cross" each other
    end

    return false if (o1 > eps && o2 > eps) || (o1 < -eps && o2 < -eps)
    return false if (o3 > eps && o4 > eps) || (o3 < -eps && o4 < -eps)


    true
  end

  # Validate a single constraint
  def validate_constraint(con, points, p1, p2)
    del = con[:del]
    x1, y1 = points[p1]
    x2, y2 = points[p2]

    found = nil
    found_adj = nil

    del[:triangles].length.times do |edg|
      e1 = del[:triangles][edg]
      e2 = del[:triangles][next_edge(edg)]

      if e1 == p1 && e2 == p2
        assert_nil found, "Constraint #{p1}->#{p2} should not be duplicate"
        found = edg
      elsif e1 == p2 && e2 == p1
        assert_nil found_adj, "Reverse constraint #{p2}->#{p1} should not be duplicate"
        found_adj = edg
      end

      # Skip edges that share an endpoint with constraint
      next if e1 == p1 || e1 == p2 || e2 == p1 || e2 == p2

      # Check no other edge intersects the constraint
      x3, y3 = points[e1]
      x4, y4 = points[e2]

      refute segments_intersect?([x1, y1], [x2, y2], [x3, y3], [x4, y4]),
        "Edge #{edg} (#{e1}->#{e2}) should not intersect constraint #{p1}->#{p2}"
    end

    assert(found || found_adj,
      "Constraint #{p1}->#{p2} should exist in triangulation")

    if found
      assert Triangles::TriVis.cdt_edge_constrained?(con, found),
        "Constraint edge #{found} should be marked constrained"
    end

    if found_adj
      assert Triangles::TriVis.cdt_edge_constrained?(con, found_adj),
        "Reverse constraint edge #{found_adj} should be marked constrained"
    end
  end

  # Validate all constraints
  def validate_all_constraints(con, points, edges)
    edges.each do |e|
      validate_constraint(con, points, e[0], e[1])
    end

    # Check that only the specified edges are constrained
    del = con[:del]
    expected_constrained = Set.new

    edges.each do |e|
      p1, p2 = e
      del[:triangles].length.times do |edg|
        e1 = del[:triangles][edg]
        e2 = del[:triangles][next_edge(edg)]
        if (e1 == p1 && e2 == p2) || (e1 == p2 && e2 == p1)
          expected_constrained << edg
        end
      end
    end

    del[:triangles].length.times do |edg|
      if expected_constrained.include?(edg)
        assert Triangles::TriVis.cdt_edge_constrained?(con, edg),
          "Edge #{edg} should be constrained"
      else
        refute Triangles::TriVis.cdt_edge_constrained?(con, edg),
          "Edge #{edg} should NOT be constrained"
      end
    end
  end

  # Exact incircle predicate — direct port of robust-predicates incircle.js
  # (MIT license, Shewchuk 1997). All macros from compile.js inlined.
  #
  # Sign convention (Shewchuk standard, positive=inside for CCW triangle):
  #   result > 0  => d is inside circumcircle of CCW(a,b,c) => Delaunay violation
  #   result <= 0 => d is outside or on circumcircle        => valid

  ICC_EPSILON    = 1.1102230246251565e-16
  ICC_SPLITTER   = 134217729.0  # 2^27 + 1
  ICC_ERRBOUND_A = (10 + 96 * ICC_EPSILON) * ICC_EPSILON

  # fast_expansion_sum_zeroelim: merge two nonoverlapping expansions
  # Both e and f must already be in non-decreasing-magnitude order (as produced
  # by scale_expansion_zeroelim).  Returns a new sorted, zero-eliminated expansion.
  def _fesze(e, elen, f, flen)
    h = []
    ei = 0; fi = 0
    enow = e[0]; fnow = f[0]
    if (fnow > enow) == (fnow > -enow)
      q = enow; ei += 1
    else
      q = fnow; fi += 1
    end
    if ei < elen && fi < flen
      enow = e[ei]; fnow = f[fi]
      if (fnow > enow) == (fnow > -enow)
        qnew = enow + q; hh = q - (qnew - enow); ei += 1
        enow = ei < elen ? e[ei] : 0.0
      else
        qnew = fnow + q; hh = q - (qnew - fnow); fi += 1
        fnow = fi < flen ? f[fi] : 0.0
      end
      q = qnew
      h << hh if hh != 0.0
      while ei < elen && fi < flen
        enow = e[ei]; fnow = f[fi]
        if (fnow > enow) == (fnow > -enow)
          x = q + enow; bvirt = x - q; hh = q - (x - bvirt) + (enow - bvirt)
          ei += 1
        else
          x = q + fnow; bvirt = x - q; hh = q - (x - bvirt) + (fnow - bvirt)
          fi += 1
        end
        q = x
        h << hh if hh != 0.0
      end
    end
    while ei < elen
      enow = e[ei]
      x = q + enow; bvirt = x - q; hh = q - (x - bvirt) + (enow - bvirt)
      q = x; ei += 1
      h << hh if hh != 0.0
    end
    while fi < flen
      fnow = f[fi]
      x = q + fnow; bvirt = x - q; hh = q - (x - bvirt) + (fnow - bvirt)
      q = x; fi += 1
      h << hh if hh != 0.0
    end
    h << q if q != 0.0 || h.empty?
    h
  end

  # scale_expansion_zeroelim: multiply expansion e by scalar b
  def _seze(e, elen, b)
    # Split b
    c = ICC_SPLITTER * b; abig = c - b; bhi = c - abig; blo = b - bhi
    h = []
    enow = e[0]
    # Two_Product_Presplit(enow, b, bhi, blo, Q, hh)
    q = enow * b
    c2 = ICC_SPLITTER * enow; abig2 = c2 - enow; ahi = c2 - abig2; alo = enow - ahi
    hh = alo * blo - (q - ahi * bhi - alo * bhi - ahi * blo)
    h << hh if hh != 0.0
    (1...elen).each do |i|
      enow = e[i]
      # Two_Product_Presplit(enow, b, bhi, blo, product1, product0)
      p1 = enow * b
      c2 = ICC_SPLITTER * enow; abig2 = c2 - enow; ahi = c2 - abig2; alo = enow - ahi
      p0 = alo * blo - (p1 - ahi * bhi - alo * bhi - ahi * blo)
      # Two_Sum(Q, product0, sum, hh)
      s = q + p0; bvirt = s - q; hh = q - (s - bvirt) + (p0 - bvirt)
      q = s
      h << hh if hh != 0.0
      # Fast_Two_Sum(product1, sum_val, Q, hh)
      qnew = p1 + q; hh = q - (qnew - p1)
      q = qnew
      h << hh if hh != 0.0
    end
    h << q if q != 0.0 || h.empty?
    h
  end

  # Cross product of 2D vectors: a*d - c*b, as exact 4-component expansion [x0,x1,x2,x3]
  # Direct expansion of $Cross_Product(a, b, c, d, D) from compile.js
  def _cross4(a, b, c, d)
    # $Two_Product(a, d, s1, s0)
    s1 = a * d
    cc = ICC_SPLITTER * a; abig = cc - a; ahi = cc - abig; alo = a - ahi
    cc = ICC_SPLITTER * d; abig = cc - d; bhi = cc - abig; blo = d - bhi
    s0 = alo * blo - (s1 - ahi * bhi - alo * bhi - ahi * blo)
    # $Two_Product(c, b, t1, t0)
    t1 = c * b
    cc = ICC_SPLITTER * c; abig = cc - c; ahi = cc - abig; alo = c - ahi
    cc = ICC_SPLITTER * b; abig = cc - b; bhi = cc - abig; blo = b - bhi
    t0 = alo * blo - (t1 - ahi * bhi - alo * bhi - ahi * blo)
    # $Two_Two_Diff(s1, s0, t1, t0, u3, D[2], D[1], D[0])
    # = $Two_One_Diff(s1, s0, t0, _j, _0, D[0])
    #     $Two_Diff(s0, t0, _i, D[0])
    _i = s0 - t0; bvirt = s0 - _i; d0 = s0 - (_i + bvirt) + (bvirt - t0)
    #     $Two_Sum(s1, _i, _j, _0)
    _j = s1 + _i; bvirt = _j - s1; _0 = s1 - (_j - bvirt) + (_i - bvirt)
    # = $Two_One_Diff(_j, _0, t1, u3, D[2], D[1])
    #     $Two_Diff(_0, t1, _i, D[1])
    _i = _0 - t1; bvirt = _0 - _i; d1 = _0 - (_i + bvirt) + (bvirt - t1)
    #     $Two_Sum(_j, _i, u3, D[2])
    u3 = _j + _i; bvirt = u3 - _j; d2 = _j - (u3 - bvirt) + (_i - bvirt)
    [d0, d1, d2, u3]
  end

  def incircle_robust(ax, ay, bx, by, cx, cy, dx, dy)
    adx = ax - dx;  ady = ay - dy
    bdx = bx - dx;  bdy = by - dy
    cdx = cx - dx;  cdy = cy - dy

    bdxcdy = bdx * cdy;  cdxbdy = cdx * bdy
    alift = adx * adx + ady * ady
    cdxady = cdx * ady;  adxcdy = adx * cdy
    blift = bdx * bdx + bdy * bdy
    adxbdy = adx * bdy;  bdxady = bdx * ady
    clift = cdx * cdx + cdy * cdy

    det = alift * (bdxcdy - cdxbdy) +
          blift * (cdxady - adxcdy) +
          clift * (adxbdy - bdxady)

    permanent = (bdxcdy.abs + cdxbdy.abs) * alift +
                (cdxady.abs + adxcdy.abs) * blift +
                (adxbdy.abs + bdxady.abs) * clift

    return det if det.abs > ICC_ERRBOUND_A * permanent

    # Adaptive exact computation
    bc = _cross4(bdx, bdy, cdx, cdy)   # bdx*cdy - cdx*bdy
    ca = _cross4(cdx, cdy, adx, ady)   # cdx*ady - adx*cdy
    ab = _cross4(adx, ady, bdx, bdy)   # adx*bdy - bdx*ady

    # adet = (adx^2 + ady^2) * bc  =>  scale bc by adx twice, by ady twice, sum
    axtbc = _seze(bc, 4, adx)
    axxbc = _seze(axtbc, axtbc.length, adx)
    aytbc = _seze(bc, 4, ady)
    ayybc = _seze(aytbc, aytbc.length, ady)
    adet  = _fesze(axxbc, axxbc.length, ayybc, ayybc.length)

    bxtca = _seze(ca, 4, bdx)
    bxxca = _seze(bxtca, bxtca.length, bdx)
    bytca = _seze(ca, 4, bdy)
    byyca = _seze(bytca, bytca.length, bdy)
    bdet  = _fesze(bxxca, bxxca.length, byyca, byyca.length)

    cxtab = _seze(ab, 4, cdx)
    cxxab = _seze(cxtab, cxtab.length, cdx)
    cytab = _seze(ab, 4, cdy)
    cyyab = _seze(cytab, cytab.length, cdy)
    cdet  = _fesze(cxxab, cxxab.length, cyyab, cyyab.length)

    abdet = _fesze(adet, adet.length, bdet, bdet.length)
    fin   = _fesze(abdet, abdet.length, cdet, cdet.length)
    fin.last
  end

  # Validate Delaunay condition for non-constrained edges
  def validate_delaunay(con, points)
    del = con[:del]
    coords = del[:coords]

    del[:triangles].length.times do |edg|
      adj = del[:half_edges][edg]
      next if Triangles::TriVis.cdt_edge_constrained?(con, edg) || adj < edg

      e1 = del[:triangles][edg]
      e2 = del[:triangles][next_edge(edg)]
      e3 = del[:triangles][next_edge(next_edge(edg))]
      a3 = del[:triangles][next_edge(next_edge(adj))]

      ax, ay = coords[e1 * 2], coords[e1 * 2 + 1]
      bx, by = coords[e2 * 2], coords[e2 * 2 + 1]
      cx, cy = coords[e3 * 2], coords[e3 * 2 + 1]
      dx, dy = coords[a3 * 2], coords[a3 * 2 + 1]

      det = Triangles::Geometry.robust_incircle(ax, ay, bx, by, cx, cy, dx, dy)

      #   det >= 0: d outside or on circumcircle = valid Delaunay
      #   det < 0:  d inside circumcircle = violation
      assert det >= 0,
        "Edge #{edg}/#{adj} violates Delaunay condition (det=#{det})"
    end
  end

  # Basic example test
  def test_example_diamond
    points = [[150, 50], [50, 200], [150, 350], [250, 200]]
    edges = [[0, 2]]

    coords = points.flatten
    del = Triangles.delaunay(coords)
    del[:coords] = coords

    con = Triangles::TriVis.cdt_init(del, edges)

    validate_halfedges(del)
    validate_vert_map(con, points)
    validate_consd_consistency(con)
    validate_all_constraints(con, points, edges)
  end

  # Generate test methods for each test file
  TEST_FILES.each do |test_file|
    test_name = "test_file_#{test_file[:name].gsub(/[^a-z0-9]/i, '_')}"

    define_method(test_name) do
      points = test_file[:points]
      edges = test_file[:edges]
      expected_error = test_file[:error]

      coords = points.flatten
      del = Triangles.delaunay(coords)
      del[:coords] = coords

      con = Triangles::TriVis.cdt_init(del)

      # Validate pre-constraint state
      validate_halfedges(del)
      validate_vert_map(con, points)

      caught_error = nil
      edges.each do |e|
        begin
          Triangles::TriVis.cdt_constrain_edge(con, e[0], e[1])
        rescue => ex
          caught_error = ex
          break unless expected_error
        end
      end

      if expected_error
        assert caught_error, "Expected error: #{expected_error}"
        assert_includes caught_error.message, expected_error,
          "Error message should contain '#{expected_error}'"
      else
        assert_nil caught_error, "Should not throw: #{caught_error&.message}"

        # Validate post-constraint state
        validate_halfedges(del)
        validate_vert_map(con, points)
        validate_consd_consistency(con)
        validate_flips_cleared(con)
        validate_all_constraints(con, points, edges)
        validate_delaunay(con, points)
      end
    end
  end
end
