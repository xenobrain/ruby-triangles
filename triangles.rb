module Triangles
  @edge_stack = Array.new(512)

  module Geometry
    # Constants for robust predicates
    ROBUST_EPSILON   = 1.1102230246251565e-16
    ROBUST_SPLITTER  = 134217729.0  # 2^27 + 1
    CCW_ERRBOUND_A   = (3 + 16 * ROBUST_EPSILON) * ROBUST_EPSILON
    CCW_ERRBOUND_B   = (2 + 12 * ROBUST_EPSILON) * ROBUST_EPSILON
    CCW_ERRBOUND_C   = (9 + 64 * ROBUST_EPSILON) * ROBUST_EPSILON * ROBUST_EPSILON
    ICC_ERRBOUND_A   = (10 + 96 * ROBUST_EPSILON) * ROBUST_EPSILON
    ICC_ERRBOUND_B   = (4 + 48 * ROBUST_EPSILON) * ROBUST_EPSILON
    ICC_ERRBOUND_C   = (44 + 576 * ROBUST_EPSILON) * ROBUST_EPSILON * ROBUST_EPSILON
    RESULT_ERRBOUND  = (3 + 8 * ROBUST_EPSILON) * ROBUST_EPSILON

    class << self
      # Signed area: positive = CCW, negative = CW, zero = collinear
      def triangle_signed_area(ax, ay, bx, by, cx, cy)
        (bx - ax) * (cy - ay) - (cx - ax) * (by - ay)
      end

      #   _seze   = scale_expansion_zeroelim  (multiply expansion by scalar)
      #   _fesze  = fast_expansion_sum_zeroelim (merge two sorted expansions)
      #   _cross4 = Two_Product difference a*d - c*b as 4-component expansion
      #   _tp_sum = Two_Product sum a*b + c*d as 4-component expansion
      #   _sq_sum = Square sum a^2 + b^2 as 4-component expansion
      def _seze(e, b)
        sp = ROBUST_SPLITTER
        cc = sp * b; abig = cc - b; bhi = cc - abig; blo = b - bhi
        h = []; enow = e[0]
        q = enow * b; cc = sp * enow; abig = cc - enow; ahi = cc - abig; alo = enow - ahi
        hh = alo * blo - (q - ahi * bhi - alo * bhi - ahi * blo); h << hh if hh != 0.0
        i = 1
        while i < e.length
          enow = e[i]; p1 = enow * b; cc = sp * enow; abig = cc - enow; ahi = cc - abig; alo = enow - ahi
          p0 = alo * blo - (p1 - ahi * bhi - alo * bhi - ahi * blo)
          s = q + p0; bvirt = s - q; hh = q - (s - bvirt) + (p0 - bvirt); q = s; h << hh if hh != 0.0
          qn = p1 + q; hh = q - (qn - p1); q = qn; h << hh if hh != 0.0
          i += 1
        end
        h << q if q != 0.0 || h.empty?; h
      end

      def _fesze(e, f)
        h = []; ei = 0; fi = 0; elen = e.length; flen = f.length
        enow = e[0]; fnow = f[0]
        if (fnow > enow) == (fnow > -enow); q = enow; ei += 1 else q = fnow; fi += 1 end
        if ei < elen && fi < flen
          enow = e[ei]; fnow = f[fi]
          if (fnow > enow) == (fnow > -enow); qn = enow + q; hh = q - (qn - enow); ei += 1 else qn = fnow + q; hh = q - (qn - fnow); fi += 1 end
          q = qn; h << hh if hh != 0.0
          while ei < elen && fi < flen
            enow = e[ei]; fnow = f[fi]
            if (fnow > enow) == (fnow > -enow); x = q + enow; bvirt = x - q; hh = q - (x - bvirt) + (enow - bvirt); ei += 1
            else; x = q + fnow; bvirt = x - q; hh = q - (x - bvirt) + (fnow - bvirt); fi += 1 end
            q = x; h << hh if hh != 0.0
          end
        end
        while ei < elen; enow = e[ei]; x = q + enow; bvirt = x - q; hh = q - (x - bvirt) + (enow - bvirt); q = x; ei += 1; h << hh if hh != 0.0 end
        while fi < flen; fnow = f[fi]; x = q + fnow; bvirt = x - q; hh = q - (x - bvirt) + (fnow - bvirt); q = x; fi += 1; h << hh if hh != 0.0 end
        h << q if q != 0.0 || h.empty?; h
      end

      def _cross4(a, b, c, d)
        sp = ROBUST_SPLITTER
        s1 = a * d; cc = sp * a; abig = cc - a; ahi = cc - abig; alo = a - ahi
        cc = sp * d; abig = cc - d; bhi = cc - abig; blo = d - bhi
        s0 = alo * blo - (s1 - ahi * bhi - alo * bhi - ahi * blo)
        t1 = c * b; cc = sp * c; abig = cc - c; ahi = cc - abig; alo = c - ahi
        cc = sp * b; abig = cc - b; bhi = cc - abig; blo = b - bhi
        t0 = alo * blo - (t1 - ahi * bhi - alo * bhi - ahi * blo)
        _i = s0 - t0; bvirt = s0 - _i; d0 = s0 - (_i + bvirt) + (bvirt - t0)
        _j = s1 + _i; bvirt = _j - s1; _z = s1 - (_j - bvirt) + (_i - bvirt)
        _i = _z - t1; bvirt = _z - _i; d1 = _z - (_i + bvirt) + (bvirt - t1)
        d3 = _j + _i; bvirt = d3 - _j; d2 = _j - (d3 - bvirt) + (_i - bvirt)
        [d0, d1, d2, d3]
      end

      def _tp_sum(a, b, c, d)
        sp = ROBUST_SPLITTER
        s1 = a * b; cc = sp * a; abig = cc - a; ahi = cc - abig; alo = a - ahi
        cc = sp * b; abig = cc - b; bhi = cc - abig; blo = b - bhi
        s0 = alo * blo - (s1 - ahi * bhi - alo * bhi - ahi * blo)
        t1 = c * d; cc = sp * c; abig = cc - c; ahi = cc - abig; alo = c - ahi
        cc = sp * d; abig = cc - d; bhi = cc - abig; blo = d - bhi
        t0 = alo * blo - (t1 - ahi * bhi - alo * bhi - ahi * blo)
        _i = s0 + t0; bvirt = _i - s0; d0 = s0 - (_i - bvirt) + (t0 - bvirt)
        _j = s1 + _i; bvirt = _j - s1; _z = s1 - (_j - bvirt) + (_i - bvirt)
        _i = _z + t1; bvirt = _i - _z; d1 = _z - (_i - bvirt) + (t1 - bvirt)
        d3 = _j + _i; bvirt = d3 - _j; d2 = _j - (d3 - bvirt) + (_i - bvirt)
        [d0, d1, d2, d3]
      end

      def _sq_sum(a, b)
        sp = ROBUST_SPLITTER
        s1 = a * a; cc = sp * a; abig = cc - a; ahi = cc - abig; alo = a - ahi
        s0 = alo * alo - (s1 - ahi * ahi - (ahi + ahi) * alo)
        t1 = b * b; cc = sp * b; abig = cc - b; ahi = cc - abig; alo = b - ahi
        t0 = alo * alo - (t1 - ahi * ahi - (ahi + ahi) * alo)
        _i = s0 + t0; bvirt = _i - s0; d0 = s0 - (_i - bvirt) + (t0 - bvirt)
        _j = s1 + _i; bvirt = _j - s1; _z = s1 - (_j - bvirt) + (_i - bvirt)
        _i = _z + t1; bvirt = _i - _z; d1 = _z - (_i - bvirt) + (t1 - bvirt)
        d3 = _j + _i; bvirt = d3 - _j; d2 = _j - (d3 - bvirt) + (_i - bvirt)
        [d0, d1, d2, d3]
      end

      def _estimate(e)
        q = e[0]; i = 1; while i < e.length; q += e[i]; i += 1; end; q
      end

      # Returns positive=CCW, negative=CW, 0=collinear (same sign as triangle_signed_area).
      # Fast path returns -det because the simple formula is negated relative to
      # the standard orient2d determinant; the adaptive path computes the standard
      # determinant directly, so no negation is needed there.
      def robust_orient2d(ax, ay, bx, by, cx, cy)
        detleft  = (ay - cy) * (bx - cx)
        detright = (ax - cx) * (by - cy)
        det = detleft - detright
        detsum = (detleft + detright).abs
        return -det if det.abs >= CCW_ERRBOUND_A * detsum
        _orient2d_adapt(ax, ay, bx, by, cx, cy, detsum)
      end

      def _orient2d_adapt(ax, ay, bx, by, cx, cy, detsum)
        sp = ROBUST_SPLITTER
        acx = ax - cx; bcx = bx - cx; acy = ay - cy; bcy = by - cy

        # Cross_Product(acx, bcx, acy, bcy) → B = acx*bcy - acy*bcx
        s1 = acx * bcy; c = sp * acx; ahi = c - (c - acx); alo = acx - ahi
        c = sp * bcy; bhi = c - (c - bcy); blo = bcy - bhi
        s0 = alo * blo - (s1 - ahi * bhi - alo * bhi - ahi * blo)
        t1 = acy * bcx; c = sp * acy; ahi = c - (c - acy); alo = acy - ahi
        c = sp * bcx; bhi = c - (c - bcx); blo = bcx - bhi
        t0 = alo * blo - (t1 - ahi * bhi - alo * bhi - ahi * blo)
        _i = s0 - t0; bvirt = s0 - _i
        b0 = s0 - (_i + bvirt) + (bvirt - t0)
        _j = s1 + _i; bvirt = _j - s1; _0 = s1 - (_j - bvirt) + (_i - bvirt)
        _i = _0 - t1; bvirt = _0 - _i
        b1 = _0 - (_i + bvirt) + (bvirt - t1)
        u3 = _j + _i; bvirt = u3 - _j
        b2 = _j - (u3 - bvirt) + (_i - bvirt)
        b3 = u3
        bb = [b0, b1, b2, b3]

        det = b0 + b1 + b2 + b3
        errbound = CCW_ERRBOUND_B * detsum
        return det if det >= errbound || -det >= errbound

        # Two_Diff_Tail for each difference
        bvirt = ax - acx; acxtail = ax - (acx + bvirt) + (bvirt - cx)
        bvirt = bx - bcx; bcxtail = bx - (bcx + bvirt) + (bvirt - cx)
        bvirt = ay - acy; acytail = ay - (acy + bvirt) + (bvirt - cy)
        bvirt = by - bcy; bcytail = by - (bcy + bvirt) + (bvirt - cy)

        return det if acxtail == 0 && acytail == 0 && bcxtail == 0 && bcytail == 0

        errbound = CCW_ERRBOUND_C * detsum + RESULT_ERRBOUND * det.abs
        det += (acx * bcytail + bcy * acxtail) - (acy * bcxtail + bcx * acytail)
        return det if det >= errbound || -det >= errbound

        u = _cross4(acxtail, bcx, acytail, bcy)
        c1 = _fesze(bb, u)

        u = _cross4(acx, bcxtail, acy, bcytail)
        c2 = _fesze(c1, u)

        u = _cross4(acxtail, bcxtail, acytail, bcytail)
        d = _fesze(c2, u)
        d.last
      end

      # Returns the incircle determinant value.
      # For CCW triangle (a,b,c): > 0 means d outside, < 0 means d inside, 0 cocircular.
      def robust_incircle(ax, ay, bx, by, cx, cy, dx, dy)
        adx = ax - dx; ady = ay - dy
        bdx = bx - dx; bdy = by - dy
        cdx = cx - dx; cdy = cy - dy

        bdxcdy = bdx * cdy; cdxbdy = cdx * bdy; alift = adx * adx + ady * ady
        cdxady = cdx * ady; adxcdy = adx * cdy; blift = bdx * bdx + bdy * bdy
        adxbdy = adx * bdy; bdxady = bdx * ady; clift = cdx * cdx + cdy * cdy

        det = alift * (bdxcdy - cdxbdy) + blift * (cdxady - adxcdy) + clift * (adxbdy - bdxady)
        permanent = (bdxcdy.abs + cdxbdy.abs) * alift +
                    (cdxady.abs + adxcdy.abs) * blift +
                    (adxbdy.abs + bdxady.abs) * clift
        errbound = ICC_ERRBOUND_A * permanent
        return det if det > errbound || -det > errbound

        _incircle_adapt(ax, ay, bx, by, cx, cy, dx, dy, permanent)
      end

      def _incircle_adapt(ax, ay, bx, by, cx, cy, dx, dy, permanent) # rubocop:disable Metrics
        adx = ax - dx; bdx = bx - dx; cdx = cx - dx
        ady = ay - dy; bdy = by - dy; cdy = cy - dy

        bc = _cross4(bdx, bdy, cdx, cdy)
        ca = _cross4(cdx, cdy, adx, ady)
        ab = _cross4(adx, ady, bdx, bdy)

        adet = _fesze(_seze(_seze(bc, adx), adx), _seze(_seze(bc, ady), ady))
        bdet = _fesze(_seze(_seze(ca, bdx), bdx), _seze(_seze(ca, bdy), bdy))
        cdet = _fesze(_seze(_seze(ab, cdx), cdx), _seze(_seze(ab, cdy), cdy))
        fin = _fesze(_fesze(adet, bdet), cdet)

        det = _estimate(fin)
        errbound = ICC_ERRBOUND_B * permanent
        return det if det >= errbound || -det >= errbound

        bvirt = ax - adx; adxtail = ax - (adx + bvirt) + (bvirt - dx)
        bvirt = ay - ady; adytail = ay - (ady + bvirt) + (bvirt - dy)
        bvirt = bx - bdx; bdxtail = bx - (bdx + bvirt) + (bvirt - dx)
        bvirt = by - bdy; bdytail = by - (bdy + bvirt) + (bvirt - dy)
        bvirt = cx - cdx; cdxtail = cx - (cdx + bvirt) + (bvirt - dx)
        bvirt = cy - cdy; cdytail = cy - (cdy + bvirt) + (bvirt - dy)

        if adxtail == 0 && bdxtail == 0 && cdxtail == 0 &&
           adytail == 0 && bdytail == 0 && cdytail == 0
          return det
        end

        errbound = ICC_ERRBOUND_C * permanent + RESULT_ERRBOUND * det.abs
        det += ((adx * adx + ady * ady) * ((bdx * cdytail + cdy * bdxtail) - (bdy * cdxtail + cdx * bdytail)) +
          2.0 * (adx * adxtail + ady * adytail) * (bdx * cdy - bdy * cdx)) +
               ((bdx * bdx + bdy * bdy) * ((cdx * adytail + ady * cdxtail) - (cdy * adxtail + adx * cdytail)) +
                 2.0 * (bdx * bdxtail + bdy * bdytail) * (cdx * ady - cdy * adx)) +
               ((cdx * cdx + cdy * cdy) * ((adx * bdytail + bdy * adxtail) - (ady * bdxtail + bdx * adytail)) +
                 2.0 * (cdx * cdxtail + cdy * cdytail) * (adx * bdy - ady * bdx))
        return det if det >= errbound || -det >= errbound

        # Full exact computation with tails
        aa = bb_sq = cc_sq = nil
        axtbc = aytbc = bxtca = bytca = cxtab = cytab = nil

        if bdxtail != 0 || bdytail != 0 || cdxtail != 0 || cdytail != 0
          aa = _sq_sum(adx, ady)
        end
        if cdxtail != 0 || cdytail != 0 || adxtail != 0 || adytail != 0
          bb_sq = _sq_sum(bdx, bdy)
        end
        if adxtail != 0 || adytail != 0 || bdxtail != 0 || bdytail != 0
          cc_sq = _sq_sum(cdx, cdy)
        end

        if adxtail != 0
          axtbc = _seze(bc, adxtail)
          t1 = _seze(axtbc, 2.0 * adx)
          t2 = _seze(_seze(cc_sq, adxtail), bdy)
          t3 = _seze(_seze(bb_sq, adxtail), -cdy)
          fin = _fesze(fin, _fesze(_fesze(t1, t2), t3))
        end
        if adytail != 0
          aytbc = _seze(bc, adytail)
          t1 = _seze(aytbc, 2.0 * ady)
          t2 = _seze(_seze(bb_sq, adytail), cdx)
          t3 = _seze(_seze(cc_sq, adytail), -bdx)
          fin = _fesze(fin, _fesze(_fesze(t1, t2), t3))
        end
        if bdxtail != 0
          bxtca = _seze(ca, bdxtail)
          t1 = _seze(bxtca, 2.0 * bdx)
          t2 = _seze(_seze(aa, bdxtail), cdy)
          t3 = _seze(_seze(cc_sq, bdxtail), -ady)
          fin = _fesze(fin, _fesze(_fesze(t1, t2), t3))
        end
        if bdytail != 0
          bytca = _seze(ca, bdytail)
          t1 = _seze(bytca, 2.0 * bdy)
          t2 = _seze(_seze(cc_sq, bdytail), adx)
          t3 = _seze(_seze(aa, bdytail), -cdx)
          fin = _fesze(fin, _fesze(_fesze(t1, t2), t3))
        end
        if cdxtail != 0
          cxtab = _seze(ab, cdxtail)
          t1 = _seze(cxtab, 2.0 * cdx)
          t2 = _seze(_seze(bb_sq, cdxtail), ady)
          t3 = _seze(_seze(aa, cdxtail), -bdy)
          fin = _fesze(fin, _fesze(_fesze(t1, t2), t3))
        end
        if cdytail != 0
          cytab = _seze(ab, cdytail)
          t1 = _seze(cytab, 2.0 * cdy)
          t2 = _seze(_seze(aa, cdytail), bdx)
          t3 = _seze(_seze(bb_sq, cdytail), -adx)
          fin = _fesze(fin, _fesze(_fesze(t1, t2), t3))
        end

        # Cross-tail terms
        if adxtail != 0 || adytail != 0
          if bdxtail != 0 || bdytail != 0 || cdxtail != 0 || cdytail != 0
            bct = _fesze(_tp_sum(bdxtail, cdy, bdx, cdytail), _tp_sum(cdxtail, -bdy, cdx, -bdytail))
            bctt = _cross4(bdxtail, bdytail, cdxtail, cdytail)
          else
            bct = [0]; bctt = [0]
          end
          if adxtail != 0
            t1 = _seze(_seze(bct, adxtail), 2.0 * adx)
            fin = _fesze(fin, _fesze(_seze(axtbc, adxtail), t1))
            t2 = _seze(_seze(bctt, adxtail), 2.0 * adx)
            t3 = _seze(_seze(bctt, adxtail), adxtail)
            t4 = _seze(_seze(bct, adxtail), adxtail)
            fin = _fesze(fin, _fesze(_fesze(t2, t3), t4))
            fin = _fesze(fin, _seze(_seze(cc_sq, adxtail), bdytail)) if bdytail != 0
            fin = _fesze(fin, _seze(_seze(bb_sq, -adxtail), cdytail)) if cdytail != 0
          end
          if adytail != 0
            t1 = _seze(_seze(bct, adytail), 2.0 * ady)
            fin = _fesze(fin, _fesze(_seze(aytbc, adytail), t1))
            t2 = _seze(_seze(bctt, adytail), 2.0 * ady)
            t3 = _seze(_seze(bctt, adytail), adytail)
            t4 = _seze(_seze(bct, adytail), adytail)
            fin = _fesze(fin, _fesze(_fesze(t2, t3), t4))
          end
        end
        if bdxtail != 0 || bdytail != 0
          if cdxtail != 0 || cdytail != 0 || adxtail != 0 || adytail != 0
            cat = _fesze(_tp_sum(cdxtail, ady, cdx, adytail), _tp_sum(adxtail, -cdy, adx, -cdytail))
            catt = _cross4(cdxtail, cdytail, adxtail, adytail)
          else
            cat = [0]; catt = [0]
          end
          if bdxtail != 0
            t1 = _seze(_seze(cat, bdxtail), 2.0 * bdx)
            fin = _fesze(fin, _fesze(_seze(bxtca, bdxtail), t1))
            t2 = _seze(_seze(catt, bdxtail), 2.0 * bdx)
            t3 = _seze(_seze(catt, bdxtail), bdxtail)
            t4 = _seze(_seze(cat, bdxtail), bdxtail)
            fin = _fesze(fin, _fesze(_fesze(t2, t3), t4))
            fin = _fesze(fin, _seze(_seze(aa, bdxtail), cdytail)) if cdytail != 0
            fin = _fesze(fin, _seze(_seze(cc_sq, -bdxtail), adytail)) if adytail != 0
          end
          if bdytail != 0
            t1 = _seze(_seze(cat, bdytail), 2.0 * bdy)
            fin = _fesze(fin, _fesze(_seze(bytca, bdytail), t1))
            t2 = _seze(_seze(catt, bdytail), 2.0 * bdy)
            t3 = _seze(_seze(catt, bdytail), bdytail)
            t4 = _seze(_seze(cat, bdytail), bdytail)
            fin = _fesze(fin, _fesze(_fesze(t2, t3), t4))
          end
        end
        if cdxtail != 0 || cdytail != 0
          if adxtail != 0 || adytail != 0 || bdxtail != 0 || bdytail != 0
            abt = _fesze(_tp_sum(adxtail, bdy, adx, bdytail), _tp_sum(bdxtail, -ady, bdx, -adytail))
            abtt = _cross4(adxtail, adytail, bdxtail, bdytail)
          else
            abt = [0]; abtt = [0]
          end
          if cdxtail != 0
            t1 = _seze(_seze(abt, cdxtail), 2.0 * cdx)
            fin = _fesze(fin, _fesze(_seze(cxtab, cdxtail), t1))
            t2 = _seze(_seze(abtt, cdxtail), 2.0 * cdx)
            t3 = _seze(_seze(abtt, cdxtail), cdxtail)
            t4 = _seze(_seze(abt, cdxtail), cdxtail)
            fin = _fesze(fin, _fesze(_fesze(t2, t3), t4))
            fin = _fesze(fin, _seze(_seze(bb_sq, cdxtail), adytail)) if adytail != 0
            fin = _fesze(fin, _seze(_seze(aa, -cdxtail), bdytail)) if bdytail != 0
          end
          if cdytail != 0
            t1 = _seze(_seze(abt, cdytail), 2.0 * cdy)
            fin = _fesze(fin, _fesze(_seze(cytab, cdytail), t1))
            t2 = _seze(_seze(abtt, cdytail), 2.0 * cdy)
            t3 = _seze(_seze(abtt, cdytail), cdytail)
            t4 = _seze(_seze(abt, cdytail), cdytail)
            fin = _fesze(fin, _fesze(_fesze(t2, t3), t4))
          end
        end

        fin.last
      end

      def segments_intersect?(p1x, p1y, p2x, p2y, p3x, p3y, p4x, p4y)
        d1 = triangle_signed_area(p3x, p3y, p4x, p4y, p1x, p1y)
        d2 = triangle_signed_area(p3x, p3y, p4x, p4y, p2x, p2y)
        d3 = triangle_signed_area(p1x, p1y, p2x, p2y, p3x, p3y)
        d4 = triangle_signed_area(p1x, p1y, p2x, p2y, p4x, p4y)

        if ((d1 > 0 && d2 < 0) || (d1 < 0 && d2 > 0)) &&
           ((d3 > 0 && d4 < 0) || (d3 < 0 && d4 > 0))
          return true
        end

        return true if d1.abs < 1e-10 && point_on_segment?(p3x, p3y, p4x, p4y, p1x, p1y)
        return true if d2.abs < 1e-10 && point_on_segment?(p3x, p3y, p4x, p4y, p2x, p2y)
        return true if d3.abs < 1e-10 && point_on_segment?(p1x, p1y, p2x, p2y, p3x, p3y)
        return true if d4.abs < 1e-10 && point_on_segment?(p1x, p1y, p2x, p2y, p4x, p4y)

        false
      end

      def point_on_segment?(ax, ay, bx, by, px, py)
        min_x = ax < bx ? ax : bx
        max_x = ax > bx ? ax : bx
        min_y = ay < by ? ay : by
        max_y = ay > by ? ay : by

        px >= min_x && px <= max_x && py >= min_y && py <= max_y
      end

      def dist_squared(ax, ay, bx, by)
        dx = ax - bx
        dy = ay - by
        dx * dx + dy * dy
      end
    end
  end

  class Node
    attr_accessor :i, :x, :y, :prev, :next, :z, :prev_z, :next_z, :steiner

    def initialize(i, x, y)
      @i = i
      @x = x
      @y = y
      @z = 0
      @steiner = false
      # @prev, @next, @prev_z, @next_z are nil by default - no need to set
    end
  end

  class << self
    # EARCUT - Polygon triangulation with holes
    # data: flat array [x0, y0, x1, y1, ...]
    # hole_indices: array of hole starts (vertex count, not coordinate count)
    # dim: coordinates per vertex (2 for 2D, 3 for 3D - only first 2 used)
    def earcut(data, hole_indices = nil, dim = 2)
      has_holes = hole_indices && hole_indices.length > 0
      outer_len = has_holes ? hole_indices[0] * dim : data.length
      outer_node = linked_list(data, 0, outer_len, dim, true)

      triangles = []
      return triangles if !outer_node || outer_node.next == outer_node.prev

      min_x = nil
      min_y = nil
      inv_size = nil

      outer_node = eliminate_holes(data, hole_indices, outer_node, dim) if has_holes

      # If the shape is not too simple, use z-order curve hash later
      if data.length > 80 * dim
        min_x = data[0]
        min_y = data[1]
        max_x = min_x
        max_y = min_y

        i = dim
        while i < outer_len
          x = data[i]
          y = data[i + 1]
          min_x = x if x < min_x
          min_y = y if y < min_y
          max_x = x if x > max_x
          max_y = y if y > max_y
          i += dim
        end

        inv_size = [max_x - min_x, max_y - min_y].max
        inv_size = inv_size != 0 ? 32767.0 / inv_size : 0
      end

      earcut_linked(outer_node, triangles, dim, min_x, min_y, inv_size, 0)
      triangles
    end

    def linked_list(data, start, end_idx, dim, clockwise)
      last = nil

      if clockwise == (signed_area(data, start, end_idx, dim) > 0)
        i = start
        while i < end_idx
          last = insert_node(i.div(dim), data[i], data[i + 1], last)
          i += dim
        end
      else
        i = end_idx - dim
        while i >= start
          last = insert_node(i.div(dim), data[i], data[i + 1], last)
          i -= dim
        end
      end

      if last && equals(last, last.next)
        remove_node(last)
        last = last.next
      end

      last
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

    def insert_node(i, x, y, last)
      p = Node.new(i, x, y)

      if last.nil?
        p.prev = p
        p.next = p
      else
        p.next = last.next
        p.prev = last
        last.next.prev = p
        last.next = p
      end

      p
    end

    def remove_node(p)
      p.next.prev = p.prev
      p.prev.next = p.next

      p.prev_z.next_z = p.next_z if p.prev_z
      p.next_z.prev_z = p.prev_z if p.next_z
    end

    def equals(p1, p2)
      p1.x == p2.x && p1.y == p2.y
    end

    # Ear-clipping with 3-pass escalation:
    #   pass 0: standard ear test (fast, handles most cases)
    #   pass 1: filter degenerate vertices and retry
    #   pass 2: cure local self-intersections by splitting polygon
    def earcut_linked(ear, triangles, dim, min_x, min_y, inv_size, pass)
      return unless ear

      index_curve(ear, min_x, min_y, inv_size) if pass == 0 && inv_size

      stop = ear

      while ear.prev != ear.next
        prev_node = ear.prev
        next_node = ear.next

        if inv_size ? ear_hashed?(ear, min_x, min_y, inv_size) : ear?(ear)
          triangles << prev_node.i << ear.i << next_node.i

          remove_node(ear)

          ear = next_node.next
          stop = next_node.next

          next
        end

        ear = next_node

        if ear == stop
          if pass == 0
            earcut_linked(filter_points(ear), triangles, dim, min_x, min_y, inv_size, 1)
          elsif pass == 1
            ear = cure_local_intersections(filter_points(ear), triangles)
            earcut_linked(ear, triangles, dim, min_x, min_y, inv_size, 2)
          elsif pass == 2
            split_earcut(ear, triangles, dim, min_x, min_y, inv_size)
          end

          break
        end
      end
    end

    def point_in_ear(p, a, c, ax, ay, bx, by, cx, cy, x0, y0, x1, y1)
      p.x >= x0 && p.x <= x1 && p.y >= y0 && p.y <= y1 && p != a && p != c &&
        !(p.x == ax && p.y == ay) && point_in_triangle(ax, ay, bx, by, cx, cy, p.x, p.y) && area(p.prev, p, p.next) >= 0
    end

    def ear?(ear)
      a = ear.prev
      c = ear.next
      return false if area(a, ear, c) >= 0

      ax, bx, cx = a.x, ear.x, c.x
      ay, by, cy = a.y, ear.y, c.y
      x0 = ax < bx ? (ax < cx ? ax : cx) : (bx < cx ? bx : cx)
      y0 = ay < by ? (ay < cy ? ay : cy) : (by < cy ? by : cy)
      x1 = ax > bx ? (ax > cx ? ax : cx) : (bx > cx ? bx : cx)
      y1 = ay > by ? (ay > cy ? ay : cy) : (by > cy ? by : cy)

      p = c.next
      while p != a
        return false if point_in_ear(p, a, c, ax, ay, bx, by, cx, cy, x0, y0, x1, y1)
        p = p.next
      end
      true
    end

    def ear_hashed?(ear, min_x, min_y, inv_size)
      a = ear.prev
      c = ear.next
      return false if area(a, ear, c) >= 0

      ax, bx, cx = a.x, ear.x, c.x
      ay, by, cy = a.y, ear.y, c.y
      x0 = ax < bx ? (ax < cx ? ax : cx) : (bx < cx ? bx : cx)
      y0 = ay < by ? (ay < cy ? ay : cy) : (by < cy ? by : cy)
      x1 = ax > bx ? (ax > cx ? ax : cx) : (bx > cx ? bx : cx)
      y1 = ay > by ? (ay > cy ? ay : cy) : (by > cy ? by : cy)

      min_z = z_order(x0, y0, min_x, min_y, inv_size)
      max_z = z_order(x1, y1, min_x, min_y, inv_size)

      p = ear.prev_z
      n = ear.next_z

      while p && p.z >= min_z && n && n.z <= max_z
        return false if point_in_ear(p, a, c, ax, ay, bx, by, cx, cy, x0, y0, x1, y1)
        p = p.prev_z
        return false if point_in_ear(n, a, c, ax, ay, bx, by, cx, cy, x0, y0, x1, y1)
        n = n.next_z
      end

      while p && p.z >= min_z
        return false if point_in_ear(p, a, c, ax, ay, bx, by, cx, cy, x0, y0, x1, y1)
        p = p.prev_z
      end

      while n && n.z <= max_z
        return false if point_in_ear(n, a, c, ax, ay, bx, by, cx, cy, x0, y0, x1, y1)
        n = n.next_z
      end

      true
    end

    def area(p, q, r)
      (q.y - p.y) * (r.x - q.x) - (q.x - p.x) * (r.y - q.y)
    end

    def point_in_triangle(ax, ay, bx, by, cx, cy, px, py)
      (cx - px) * (ay - py) >= (ax - px) * (cy - py) &&
        (ax - px) * (by - py) >= (bx - px) * (ay - py) &&
        (bx - px) * (cy - py) >= (cx - px) * (by - py)
    end

    def filter_points(start, ending = nil)
      return start unless start
      ending ||= start

      p = start
      again = nil

      while true
        again = false

        if !p.steiner && (equals(p, p.next) || area(p.prev, p, p.next) == 0)
          remove_node(p)
          p = ending = p.prev
          break if p == p.next
          again = true
        else
          p = p.next
        end

        break unless again || p != ending
      end

      ending
    end

    def cure_local_intersections(start, triangles)
      p = start
      while true
        a = p.prev
        b = p.next.next

        if !equals(a, b) && intersects(a, p, p.next, b) && locally_inside(a, b) && locally_inside(b, a)
          triangles << a.i << p.i << b.i

          remove_node(p)
          remove_node(p.next)

          p = start = b
        end
        p = p.next
        break if p == start
      end

      filter_points(p)
    end

    def split_earcut(start, triangles, dim, min_x, min_y, inv_size)
      a = start
      while true
        b = a.next.next
        while b != a.prev
          if a.i != b.i && valid_diagonal?(a, b)
            c = split_polygon(a, b)

            a = filter_points(a, a.next)
            c = filter_points(c, c.next)

            earcut_linked(a, triangles, dim, min_x, min_y, inv_size, 0)
            earcut_linked(c, triangles, dim, min_x, min_y, inv_size, 0)
            return
          end
          b = b.next
        end
        a = a.next
        break if a == start
      end
    end

    def eliminate_holes(data, hole_indices, outer_node, dim)
      queue = []
      len = hole_indices.length

      i = 0
      while i < len
        start = hole_indices[i] * dim
        end_idx = i < len - 1 ? hole_indices[i + 1] * dim : data.length
        list = linked_list(data, start, end_idx, dim, false)
        list.steiner = true if list == list.next
        queue << get_leftmost(list)
        i += 1
      end

      queue.sort! { |a, b| compare_x_y_slope(a, b) }

      i = 0
      while i < queue.length
        outer_node = eliminate_hole(queue[i], outer_node)
        i += 1
      end

      outer_node
    end

    def compare_x_y_slope(a, b)
      result = a.x - b.x
      if result == 0
        result = a.y - b.y
        if result == 0
          a_slope = (a.next.y - a.y).to_f / (a.next.x - a.x)
          b_slope = (b.next.y - b.y).to_f / (b.next.x - b.x)
          result = a_slope - b_slope
          return 0 if result.nan?
        end
      end
      result <=> 0
    end

    def eliminate_hole(hole, outer_node)
      bridge = find_hole_bridge(hole, outer_node)
      return outer_node unless bridge

      bridge_reverse = split_polygon(bridge, hole)

      filter_points(bridge_reverse, bridge_reverse.next)
      filter_points(bridge, bridge.next)
    end

    def find_hole_bridge(hole, outer_node)
      p = outer_node
      hx = hole.x
      hy = hole.y
      qx = -Float::INFINITY
      m = nil

      return p if equals(hole, p)
      while true
        return p.next if equals(hole, p.next)
        if hy <= p.y && hy >= p.next.y && p.next.y != p.y
          x = p.x + (hy - p.y) * (p.next.x - p.x).to_f / (p.next.y - p.y)
          if x <= hx && x > qx
            qx = x
            m = p.x < p.next.x ? p : p.next
            return m if x == hx
          end
        end
        p = p.next
        break if p == outer_node
      end

      return nil unless m

      stop = m
      mx = m.x
      my = m.y
      tan_min = Float::INFINITY

      p = m

      while true
        if hx >= p.x && p.x >= mx && hx != p.x &&
           point_in_triangle(hy < my ? hx : qx, hy, mx, my, hy < my ? qx : hx, hy, p.x, p.y)

          tan = (hy - p.y).abs.to_f / (hx - p.x)

          if locally_inside(p, hole) &&
             (tan < tan_min || (tan == tan_min && (p.x > m.x || (p.x == m.x && sector_contains_sector(m, p)))))
            m = p
            tan_min = tan
          end
        end

        p = p.next
        break if p == stop
      end

      m
    end

    def sector_contains_sector(m, p)
      area(m.prev, m, p.prev) < 0 && area(p.next, m, m.next) < 0
    end

    # Find the leftmost node of a polygon ring
    def get_leftmost(start)
      p = start
      leftmost = start
      while true
        leftmost = p if p.x < leftmost.x || (p.x == leftmost.x && p.y < leftmost.y)
        p = p.next
        break if p == start
      end
      leftmost
    end

    # Check if a diagonal between two polygon nodes is valid
    def valid_diagonal?(a, b)
      a.next.i != b.i && a.prev.i != b.i && !intersects_polygon(a, b) &&
        ((locally_inside(a, b) && locally_inside(b, a) && middle_inside(a, b) &&
          (area(a.prev, a, b.prev) != 0 || area(a, b.prev, b) != 0)) ||
          (equals(a, b) && area(a.prev, a, a.next) > 0 && area(b.prev, b, b.next) > 0))
    end

    def intersects(p1, q1, p2, q2)
      Geometry.segments_intersect?(p1.x, p1.y, q1.x, q1.y, p2.x, p2.y, q2.x, q2.y)
    end

    def on_segment(p, q, r)
      Geometry.point_on_segment?(p.x, p.y, r.x, r.y, q.x, q.y)
    end

    def intersects_polygon(a, b)
      p = a
      while true
        if p.i != a.i && p.next.i != a.i && p.i != b.i && p.next.i != b.i && intersects(p, p.next, a, b)
          return true
        end
        p = p.next
        break if p == a
      end
      false
    end

    def locally_inside(a, b)
      if area(a.prev, a, a.next) < 0
        area(a, b, a.next) >= 0 && area(a, a.prev, b) >= 0
      else
        area(a, b, a.prev) < 0 || area(a, a.next, b) < 0
      end
    end

    def middle_inside(a, b)
      p = a
      inside = false
      px = (a.x + b.x) / 2.0
      py = (a.y + b.y) / 2.0
      while true
        if ((p.y > py) != (p.next.y > py)) && p.next.y != p.y &&
           (px < (p.next.x - p.x) * (py - p.y) / (p.next.y - p.y) + p.x)
          inside = !inside
        end
        p = p.next
        break if p == a
      end
      inside
    end

    def split_polygon(a, b)
      a2 = Node.new(a.i, a.x, a.y)
      b2 = Node.new(b.i, b.x, b.y)
      an = a.next
      bp = b.prev

      a.next = b
      b.prev = a

      a2.next = an
      an.prev = a2

      b2.next = a2
      a2.prev = b2

      bp.next = b2
      b2.prev = bp

      b2
    end

    def index_curve(start, min_x, min_y, inv_size)
      p = start
      while true
        p.z = z_order(p.x, p.y, min_x, min_y, inv_size) if p.z == 0
        p.prev_z = p.prev
        p.next_z = p.next
        p = p.next
        break if p == start
      end

      p.prev_z.next_z = nil
      p.prev_z = nil

      sort_linked(p)
    end

    def sort_linked(list)
      in_size = 1

      while true
        p = list
        list = nil
        tail = nil
        num_merges = 0

        while p
          num_merges += 1
          q = p
          p_size = 0

          i = 0
          while i < in_size
            p_size += 1
            q = q.next_z
            break unless q
            i += 1
          end

          q_size = in_size

          while p_size > 0 || (q_size > 0 && q)
            if p_size != 0 && (q_size == 0 || !q || p.z <= q.z)
              e = p
              p = p.next_z
              p_size -= 1
            else
              e = q
              q = q.next_z
              q_size -= 1
            end

            if tail
              tail.next_z = e
            else
              list = e
            end

            e.prev_z = tail
            tail = e
          end

          p = q
        end

        tail.next_z = nil
        in_size *= 2

        break if num_merges <= 1
      end

      list
    end

    def z_order(x, y, min_x, min_y, inv_size)
      x = ((x - min_x) * inv_size).to_i
      y = ((y - min_y) * inv_size).to_i

      x = (x | (x << 8)) & 0x00FF00FF
      x = (x | (x << 4)) & 0x0F0F0F0F
      x = (x | (x << 2)) & 0x33333333
      x = (x | (x << 1)) & 0x55555555

      y = (y | (y << 8)) & 0x00FF00FF
      y = (y | (y << 4)) & 0x0F0F0F0F
      y = (y | (y << 2)) & 0x33333333
      y = (y | (y << 1)) & 0x55555555

      x | (y << 1)
    end

    # QUICKHULL - Convex hull computation
    # coords: flat array [x0, y0, x1, y1, ...]
    # tol: tolerance for collinear points (default: 0.0)
    # Returns: indices forming convex hull in CCW order
    def convex_hull(coords, tol = 0.0)
      n = coords.length.div(2)
      return [] if n == 0
      return [0] if n == 1
      return [0, 1] if n == 2

      start_idx = 0
      end_idx = 0
      min_x = coords[0]
      min_y = coords[1]
      max_x = min_x
      max_y = min_y

      i = 1
      while i < n
        x = coords[2 * i]
        y = coords[2 * i + 1]

        if x < min_x || (x == min_x && y < min_y)
          min_x = x
          min_y = y
          start_idx = i
        end

        if x > max_x || (x == max_x && y > max_y)
          max_x = x
          max_y = y
          end_idx = i
        end

        i += 1
      end

      return [start_idx] if start_idx == end_idx

      work_indices = []
      i = 0
      while i < n
        unless i == start_idx || i == end_idx
          work_indices << i
        end
        i += 1
      end

      result = [start_idx]

      # Partition and reduce on the upper hull (from start to end)
      ax = coords[2 * start_idx]
      ay = coords[2 * start_idx + 1]
      bx = coords[2 * end_idx]
      by = coords[2 * end_idx + 1]

      upper_points = qhull_partition(coords, work_indices, ax, ay, bx, by, tol)
      qhull_reduce(coords, upper_points, ax, ay, bx, by, tol, result)

      result << end_idx

      lower_points = qhull_partition(coords, work_indices, bx, by, ax, ay, tol)
      qhull_reduce(coords, lower_points, bx, by, ax, ay, tol, result)

      result
    end

    def qhull_partition(coords, indices, ax, ay, bx, by, tol)
      return [] if indices.empty?

      max_dist = 0.0
      pivot_idx = -1

      dx = bx - ax
      dy = by - ay
      line_length = Math.sqrt(dx * dx + dy * dy)
      value_tol = tol * line_length

      left_points = []

      ii = 0
      while ii < indices.length
        idx = indices[ii]
        px = coords[2 * idx]
        py = coords[2 * idx + 1]

        cross = (px - ax) * dy - (py - ay) * dx

        if cross > value_tol
          left_points << idx

          if cross > max_dist
            max_dist = cross
            pivot_idx = idx
          end
        end
        ii += 1
      end

      if pivot_idx >= 0 && left_points[0] != pivot_idx
        left_points.delete(pivot_idx)
        left_points.unshift(pivot_idx)
      end

      left_points
    end

    def qhull_reduce(coords, indices, ax, ay, bx, by, tol, result)
      return if indices.empty?

      pivot_idx = indices[0]
      px = coords[2 * pivot_idx]
      py = coords[2 * pivot_idx + 1]

      if indices.length == 1
        result << pivot_idx
        return
      end

      remaining = indices[1..-1]

      left_points = qhull_partition(coords, remaining, ax, ay, px, py, tol)
      qhull_reduce(coords, left_points, ax, ay, px, py, tol, result)

      result << pivot_idx

      right_points = qhull_partition(coords, remaining, px, py, bx, by, tol)
      qhull_reduce(coords, right_points, px, py, bx, by, tol, result)
    end

    # Delaunay triangulation.
    # coords: flat [x0,y0,x1,y1,...].
    # Returns {triangles:, half_edges:, hull:} (half-edge data structure).
    def delaunay(coords)
      n = coords.length >> 1

      return { hull: [], triangles: [], half_edges: [] } if n == 0
      return { hull: [0], triangles: [], half_edges: [] } if n == 1
      if n == 2
        if coords[0] < coords[2] || (coords[0] == coords[2] && coords[1] < coords[3])
          return { hull: [0, 1], triangles: [], half_edges: [] }
        else
          return { hull: [1, 0], triangles: [], half_edges: [] }
        end
      end

      max_triangles = 2 * n - 5
      triangles = Array.new(max_triangles * 3)
      half_edges = Array.new(max_triangles * 3)

      link = proc do |a, b|
        half_edges[a] = b
        half_edges[b] = a if b != -1
      end

      hash_size = Math.sqrt(n).ceil
      hull_prev = Array.new(n)
      hull_next = Array.new(n)
      hull_tri = Array.new(n, 0)
      hull_hash = Array.new(hash_size)

      ids = Array.new(n)
      dists = Array.new(n)

      min_x = Float::INFINITY
      min_y = Float::INFINITY
      max_x = -Float::INFINITY
      max_y = -Float::INFINITY

      i = 0
      while i < n
        x = coords[2 * i]
        y = coords[2 * i + 1]

        min_x = x if x < min_x
        min_y = y if y < min_y
        max_x = x if x > max_x
        max_y = y if y > max_y
        ids[i] = i
        i += 1
      end

      cx = (min_x + max_x) * 0.5
      cy = (min_y + max_y) * 0.5

      hash_key = proc { |kx, ky| (pseudo_angle(kx - cx, ky - cy) * hash_size).floor % hash_size }

      i0 = i1 = i2 = 0
      min_dist = Float::INFINITY

      i = 0
      while i < n
        d = Geometry.dist_squared(cx, cy, coords[2 * i], coords[2 * i + 1])
        if d < min_dist
          i0 = i
          min_dist = d
        end

        i += 1
      end

      i0x = coords[2 * i0]
      i0y = coords[2 * i0 + 1]
      min_dist = Float::INFINITY

      i = 0
      while i < n
        if i == i0
          i += 1
          next
        end

        d = Geometry.dist_squared(i0x, i0y, coords[2 * i], coords[2 * i + 1])
        if d < min_dist && d > 0
          i1 = i
          min_dist = d
        end

        i += 1
      end

      i1x = coords[2 * i1]
      i1y = coords[2 * i1 + 1]
      min_radius = Float::INFINITY

      i = 0
      while i < n
        if i == i0 || i == i1
          i += 1
          next
        end
        r = circumradius(i0x, i0y, i1x, i1y, coords[2 * i], coords[2 * i + 1])
        if r < min_radius
          i2 = i
          min_radius = r
        end

        i += 1
      end

      i2x = coords[2 * i2]
      i2y = coords[2 * i2 + 1]

      if min_radius == Float::INFINITY
        i = 0
        while i < n
          dists[i] = (coords[2 * i] - coords[0]).nonzero? || (coords[2 * i + 1] - coords[1])
          i += 1
        end

        ids = ids.sort_by { |id| dists[id] }

        hull = Array.new(n)
        d0 = -Float::INFINITY
        j = 0
        i = 0
        while i < n
          id = ids[i]
          d = dists[id]
          if d > d0
            hull[j] = id
            d0 = d
            j += 1
          end
          i += 1
        end

        hull = hull.slice(0, j)
        triangles = Array.new(0)
        half_edges = Array.new(0)
        return { hull: hull, triangles: triangles, half_edges: half_edges }
      end

      if Geometry.robust_orient2d(i0x, i0y, i1x, i1y, i2x, i2y) > 0
        i = i1
        x = i1x
        y = i1y
        i1 = i2
        i1x = i2x
        i1y = i2y
        i2 = i
        i2x = x
        i2y = y
      end

      center = circumcenter(i0x, i0y, i1x, i1y, i2x, i2y)
      cx, cy = center

      i = 0
      while i < n
        dists[i] = Geometry.dist_squared(coords[2 * i], coords[2 * i + 1], cx, cy)

        i += 1
      end

      ids = ids.sort_by { |id| dists[id] }

      hull_start = i0
      hull_size = 3

      hull_next[i0] = hull_prev[i2] = i1
      hull_next[i1] = hull_prev[i0] = i2
      hull_next[i2] = hull_prev[i1] = i0
      hull_tri[i0] = 0
      hull_tri[i1] = 1
      hull_tri[i2] = 2

      i = 0
      while i < hash_size
        hull_hash[i] = -1
        i += 1
      end
      hull_hash[hash_key[i0x, i0y]] = i0
      hull_hash[hash_key[i1x, i1y]] = i1
      hull_hash[hash_key[i2x, i2y]] = i2

      triangles_len = 0

      add_triangle = proc do |ti0, ti1, ti2, a, b, c|
        t = triangles_len

        triangles[t] = ti0
        triangles[t + 1] = ti1
        triangles[t + 2] = ti2

        link[t, a]
        link[t + 1, b]
        link[t + 2, c]

        triangles_len += 3
        t
      end

      legalize = proc do |a|
        stack_idx = ar = 0

        while true
          b = half_edges[a]

          a0 = a - a % 3
          ar = a0 + (a + 2) % 3

          if b == -1
            break if stack_idx == 0
            stack_idx -= 1
            a = @edge_stack[stack_idx]
            next
          end

          b0 = b - b % 3
          al = a0 + (a + 1) % 3
          bl = b0 + (b + 2) % 3

          p0 = triangles[ar]
          pr = triangles[a]
          pl = triangles[al]
          p1 = triangles[bl]

          illegal = in_circle(
            coords[2 * p0], coords[2 * p0 + 1],
            coords[2 * pr], coords[2 * pr + 1],
            coords[2 * pl], coords[2 * pl + 1],
            coords[2 * p1], coords[2 * p1 + 1]
          )

          if illegal
            triangles[a] = p1
            triangles[b] = p0

            hbl = half_edges[bl]

            if hbl == -1
              hull_e = hull_start
              while true
                if hull_tri[hull_e] == bl
                  hull_tri[hull_e] = a
                  break
                end
                hull_e = hull_prev[hull_e]
                break if hull_e == hull_start
              end
            end

            link[a, hbl]
            link[b, half_edges[ar]]
            link[ar, bl]

            br = b0 + (b + 1) % 3

            if stack_idx < @edge_stack.length
              @edge_stack[stack_idx] = br
              stack_idx += 1
            end
          else
            break if stack_idx == 0
            stack_idx -= 1
            a = @edge_stack[stack_idx]
          end
        end

        ar
      end

      add_triangle[i0, i1, i2, -1, -1, -1]
      xp = yp = 0.0
      k = 0

      while k < n
        i = ids[k]
        x = coords[2 * i]
        y = coords[2 * i + 1]

        if k > 0 && (x - xp).abs <= Float::EPSILON && (y - yp).abs <= Float::EPSILON
          k += 1
          next
        end
        xp = x
        yp = y

        if i == i0 || i == i1 || i == i2
          k += 1
          next
        end

        start = 0
        j = 0
        key = hash_key[x, y]
        while j < hash_size
          start = hull_hash[(key + j) % hash_size]
          break if start != -1 && start != hull_next[start]
          j += 1
        end

        start = hull_prev[start]
        e = start

        while Geometry.robust_orient2d(x, y, coords[2 * e], coords[2 * e + 1], coords[2 * hull_next[e]], coords[2 * hull_next[e] + 1]) <= 0
          e = hull_next[e]
          if e == start
            e = -1
            break
          end
        end

        if e == -1
          k += 1
          next
        end

        t = add_triangle[e, i, hull_next[e], -1, -1, hull_tri[e]]

        hull_tri[i] = legalize[t + 2]
        hull_tri[e] = t
        hull_size += 1

        next_hull = hull_next[e]
        q = hull_next[next_hull]
        while Geometry.robust_orient2d(x, y, coords[2 * next_hull], coords[2 * next_hull + 1], coords[2 * q], coords[2 * q + 1]) > 0
          t = add_triangle[next_hull, i, q, hull_tri[i], -1, hull_tri[next_hull]]
          hull_tri[i] = legalize[t + 2]
          hull_next[next_hull] = next_hull
          hull_size -= 1
          next_hull = q
          q = hull_next[next_hull]
        end

        if e == start
          q = hull_prev[e]
          while Geometry.robust_orient2d(x, y, coords[2 * q], coords[2 * q + 1], coords[2 * e], coords[2 * e + 1]) > 0
            t = add_triangle[q, i, e, -1, hull_tri[e], hull_tri[q]]
            legalize[t + 2]
            hull_tri[q] = t
            hull_next[e] = e
            hull_size -= 1
            e = q
            q = hull_prev[e]
          end
        end

        hull_start = hull_prev[i] = e
        hull_next[e] = hull_prev[next_hull] = i
        hull_next[i] = next_hull

        hull_hash[hash_key[x, y]] = i
        hull_hash[hash_key[coords[2 * e], coords[2 * e + 1]]] = e

        k += 1
      end

      hull = Array.new(hull_size)
      e = hull_start
      i = 0
      while i < hull_size
        hull[i] = e
        e = hull_next[e]
        i += 1
      end

      triangles = triangles.slice(0, triangles_len)
      half_edges = half_edges.slice(0, triangles_len)

      { hull: hull, triangles: triangles, half_edges: half_edges }
    end

    def pseudo_angle(dx, dy)
      p = dx / (dx.abs + dy.abs)
      (dy > 0 ? 3 - p : 1 + p) / 4
    end

    def in_circle(ax, ay, bx, by, cx, cy, px, py)
      dx = ax - px
      dy = ay - py
      ex = bx - px
      ey = by - py
      fx = cx - px
      fy = cy - py

      ap = dx * dx + dy * dy
      bp = ex * ex + ey * ey
      cp = fx * fx + fy * fy

      dx * (ey * cp - bp * fy) -
        dy * (ex * cp - bp * fx) +
        ap * (ex * fy - ey * fx) < 0
    end

    def circumdelta(ax, ay, bx, by, cx, cy)
      dx = bx - ax
      dy = by - ay
      ex = cx - ax
      ey = cy - ay
      bl = dx * dx + dy * dy
      cl = ex * ex + ey * ey
      d = 0.5 / (dx * ey - dy * ex)
      [(ey * bl - dy * cl) * d, (dx * cl - ex * bl) * d]
    end

    def circumradius(ax, ay, bx, by, cx, cy)
      x, y = circumdelta(ax, ay, bx, by, cx, cy)
      x * x + y * y
    end

    def circumcenter(ax, ay, bx, by, cx, cy)
      x, y = circumdelta(ax, ay, bx, by, cx, cy)
      [ax + x, ay + y]
    end

    # Polygon utility methods — all take flat coordinate arrays [x0, y0, x1, y1, ...]

    # Signed polygon area via shoelace formula. Positive = CCW, negative = CW.
    def polygon_area(coords)
      n = coords.length.div(2)
      return 0.0 if n < 3

      sum = 0.0
      i = 0
      while i < n
        j = (i + 1) % n
        x1 = coords[i * 2];     y1 = coords[i * 2 + 1]
        x2 = coords[j * 2];     y2 = coords[j * 2 + 1]
        sum += (x1 * y2) - (x2 * y1)
        i += 1
      end

      sum * 0.5
    end

    # Ensure polygon winds clockwise (negative area). Returns flat array.
    def ensure_clockwise(coords)
      polygon_area(coords) > 0 ? reverse_polygon(coords) : coords
    end

    # Ensure polygon winds counter-clockwise (positive area). Returns flat array.
    def ensure_ccw(coords)
      polygon_area(coords) < 0 ? reverse_polygon(coords) : coords
    end

    # Centroid of a simple polygon. Returns [cx, cy].
    def centroid(coords)
      area = polygon_area(coords)
      return [0.0, 0.0] if area == 0.0

      cx = 0.0
      cy = 0.0
      n = coords.length.div(2)
      i = 0
      while i < n
        j = (i + 1) % n
        x1 = coords[i * 2];     y1 = coords[i * 2 + 1]
        x2 = coords[j * 2];     y2 = coords[j * 2 + 1]
        cross = (x1 * y2) - (x2 * y1)
        cx += (x1 + x2) * cross
        cy += (y1 + y2) * cross
        i += 1
      end

      factor = 1.0 / (6.0 * area)
      [cx * factor, cy * factor]
    end

    # Ray-casting point-in-polygon test. coords is flat [x0,y0,x1,y1,...].
    def point_in_poly?(x, y, coords)
      inside = false
      n = coords.length.div(2)
      j = n - 1
      i = 0
      while i < n
        xi = coords[i * 2];     yi = coords[i * 2 + 1]
        xj = coords[j * 2];     yj = coords[j * 2 + 1]

        if ((yi > y) != (yj > y)) && (x < (xj - xi) * (y - yi) / (yj - yi) + xi)
          inside = !inside
        end

        j = i
        i += 1
      end
      inside
    end

    # Axis-aligned bounding box. Returns [min_x, min_y, max_x, max_y].
    def bounding_box(coords)
      return [0.0, 0.0, 0.0, 0.0] if coords.empty?

      min_x = Float::INFINITY
      min_y = Float::INFINITY
      max_x = -Float::INFINITY
      max_y = -Float::INFINITY

      i = 0
      while i < coords.length
        x = coords[i]; y = coords[i + 1]
        min_x = x if x < min_x
        min_y = y if y < min_y
        max_x = x if x > max_x
        max_y = y if y > max_y
        i += 2
      end

      [min_x, min_y, max_x, max_y]
    end

    # Test whether a polygon is convex. coords are flat [x0,y0,...].
    def convex?(coords)
      n = coords.length.div(2)
      return false if n < 3

      sign = nil
      i = 0
      while i < n
        j = (i + 1) % n
        k = (i + 2) % n
        x1 = coords[i * 2];     y1 = coords[i * 2 + 1]
        x2 = coords[j * 2];     y2 = coords[j * 2 + 1]
        x3 = coords[k * 2];     y3 = coords[k * 2 + 1]

        cross = (x2 - x1) * (y3 - y2) - (y2 - y1) * (x3 - x2)
        next i += 1 if cross.abs < 1e-10

        current_sign = cross > 0
        if sign.nil?
          sign = current_sign
        elsif sign != current_sign
          return false
        end
        i += 1
      end
      true
    end

    # Decompose non-convex polygon into convex parts (Bayazit algorithm).
    # coords: flat [x0,y0,...]. Returns array of flat arrays.
    def decompose_to_convex_hulls(coords)
      n = coords.length.div(2)
      return [coords] if n < 3
      return [coords] if convex?(coords)

      input_area = polygon_area(coords)
      input_is_ccw = input_area > 0

      vertices = coords.dup
      vertices = reverse_polygon(vertices) if polygon_area(vertices) < 0

      result = []
      Shatter.bayazit_decompose(vertices, result)

      unless input_is_ccw
        i = 0
        while i < result.length
          result[i] = reverse_polygon(result[i])
          i += 1
        end
      end
      result
    end

    # Shatter polygon into fragments using Voronoi-based fracturing.
    # coords: flat [x0,y0,...]. focus: [fx,fy]. Returns array of flat arrays.
    def shatter_polygon(coords, focus:, cell_size:, jitter: 0.45, min_area: 40.0, focus_radius: nil,
                        focus_extra_sites: 0, outer_site_chance: 1.0)
      Shatter.shatter_polygon(coords, focus: focus, cell_size: cell_size, jitter: jitter,
                              min_area: min_area, focus_radius: focus_radius,
                              focus_extra_sites: focus_extra_sites, outer_site_chance: outer_site_chance)
    end

    # Reverse vertex order in a flat coordinate array, preserving xy pairs.
    def reverse_polygon(coords)
      n = coords.length.div(2)
      result = Array.new(coords.length)
      i = 0
      while i < n
        j = n - 1 - i
        result[i * 2]     = coords[j * 2]
        result[i * 2 + 1] = coords[j * 2 + 1]
        i += 1
      end
      result
    end
  end
end

module Triangles
  # Polygon fracturing and geometry utilities
  module Shatter
    class << self

      # Generate Worley noise sites within a polygon. coords: flat [x0,y0,...].
      # Sites are returned as [x,y] pairs (single points, not polygons).
      def generate_worley_sites(coords, cell_size:, jitter:, focus: nil, seed: nil, focus_radius: nil,
                                focus_extra_sites: 0, outer_site_chance: 1.0)
        jitter = [[jitter, 0.0].max, 0.5].min
        border = 0.5 - jitter
        seed = rand(0x1_0000_0000) if seed.nil?

        focus_radius_sq = focus_radius ? focus_radius * focus_radius : nil

        min_x, min_y, max_x, max_y = Triangles.bounding_box(coords)
        width = max_x - min_x
        height = max_y - min_y

        cols = [(width / cell_size).floor + 1, 1].max
        rows = [(height / cell_size).floor + 1, 1].max

        center_x = (min_x + max_x) * 0.5
        center_y = (min_y + max_y) * 0.5

        all_sites = []
        inside_flags = []
        inside_sites = []

        if focus && Triangles.point_in_poly?(focus[0], focus[1], coords)
          all_sites << focus
          inside_flags << true
          inside_sites << focus
        end

        row = 0
        while row < rows
          col = 0
          while col < cols
            h = ((col * 1_640_531_513) ^ (row * 2_654_435_789)) + seed
            h &= 0xFFFF_FFFF
            span = 1.0 - border * 2.0
            fx = border + span * ((h & 0xFFFF) / 65_535.0)
            fy = border + span * (((h >> 16) & 0xFFFF) / 65_535.0)
            x = center_x + cell_size * (col + fx - cols * 0.5)
            y = center_y + cell_size * (row + fy - rows * 0.5)

            skip_site = false
            if focus_radius_sq && focus
              dx = x - focus[0]
              dy = y - focus[1]
              skip_site = (dx * dx + dy * dy) > focus_radius_sq && rand > outer_site_chance
            end

            unless skip_site
              inside = Triangles.point_in_poly?(x, y, coords)

              all_sites << [x, y]
              inside_flags << inside
              inside_sites << [x, y] if inside
            end

            col += 1
          end
          row += 1
        end

        if focus && focus_radius && focus_extra_sites.to_i > 0
          i = 0
          num_extra = focus_extra_sites.to_i
          while i < num_extra
            angle = rand * Math::PI * 2.0
            radius = Math.sqrt(rand) * focus_radius
            x = focus[0] + Math.cos(angle) * radius
            y = focus[1] + Math.sin(angle) * radius
            inside = Triangles.point_in_poly?(x, y, coords)

            all_sites << [x, y]
            inside_flags << inside
            inside_sites << [x, y] if inside
            i += 1
          end
        end

        inside_sites = [Triangles.centroid(coords)] if inside_sites.empty?

        {
          inside: inside_sites,
          all: all_sites,
          inside_flags: inside_flags
        }
      end

      # Clip a polygon by a half-plane. coords: flat [x0,y0,...].
      # Keeps the side where dot(point, [nx,ny]) <= threshold.
      # Returns flat array.
      def clip_polygon_half_plane(coords, nx, ny, threshold)
        n = coords.length.div(2)
        return [] if n < 3

        output = []
        prev_x = coords[(n - 1) * 2]
        prev_y = coords[(n - 1) * 2 + 1]
        prev_dist = (prev_x * nx + prev_y * ny) - threshold

        i = 0
        while i < n
          curr_x = coords[i * 2]
          curr_y = coords[i * 2 + 1]
          curr_dist = (curr_x * nx + curr_y * ny) - threshold

          if curr_dist <= 0.0
            if prev_dist > 0.0
              t = prev_dist / (prev_dist - curr_dist)
              output << prev_x + (curr_x - prev_x) * t
              output << prev_y + (curr_y - prev_y) * t
            end
            output << curr_x << curr_y
          elsif prev_dist <= 0.0
            t = prev_dist / (prev_dist - curr_dist)
            output << prev_x + (curr_x - prev_x) * t
            output << prev_y + (curr_y - prev_y) * t
          end

          prev_x = curr_x
          prev_y = curr_y
          prev_dist = curr_dist
          i += 1
        end

        output
      end

      # Voronoi-based polygon fracturing. coords: flat [x0,y0,...].
      # focus: [fx,fy]. Returns array of flat arrays.
      def shatter_polygon(coords, focus:, cell_size:, jitter: 0.45, min_area: 40.0, focus_radius: nil,
                          focus_extra_sites: 0, outer_site_chance: 1.0)
        site_sets = generate_worley_sites(
          coords,
          cell_size: cell_size,
          jitter: jitter,
          focus: focus,
          focus_radius: focus_radius,
          focus_extra_sites: focus_extra_sites,
          outer_site_chance: outer_site_chance
        )
        inside_sites = site_sets[:inside]
        all_sites = site_sets[:all]
        inside_flags = site_sets[:inside_flags]
        shards = []

        si = 0
        while si < inside_sites.length
          site = inside_sites[si]
          cell = coords

          ai = 0
          while ai < all_sites.length
            other = all_sites[ai]
            unless other.equal?(site) || !inside_flags[ai]
              nx = other[0] - site[0]
              ny = other[1] - site[1]
              threshold = ((site[0] + other[0]) * 0.5 * nx) +
                          ((site[1] + other[1]) * 0.5 * ny)
              cell = clip_polygon_half_plane(cell, nx, ny, threshold)
              break if cell.length < 6
            end
            ai += 1
          end

          area = Triangles.polygon_area(cell).abs
          shards << Triangles.ensure_ccw(cell) if cell.length >= 6 && area >= min_area
          si += 1
        end

        shards
      end

      # Bayazit convex decomposition. coords: flat [x0,y0,...].
      # Appends convex flat arrays to list.
      def bayazit_decompose(coords, list)
        n = coords.length.div(2)
        lower_int = [0.0, 0.0]
        upper_int = [0.0, 0.0]
        lower_index = 0
        upper_index = 0

        i = 0
        while i < n
          if bayazit_reflex?(i, coords)
            lower_dist = Float::INFINITY
            upper_dist = Float::INFINITY

            j = 0
            while j < n
              vi_1x, vi_1y = bayazit_at(i - 1, coords)
              vix,   viy   = bayazit_at(i, coords)
              vjx,   vjy   = bayazit_at(j, coords)
              vj_1x, vj_1y = bayazit_at(j - 1, coords)

              if Geometry.triangle_signed_area(vi_1x, vi_1y, vix, viy, vjx, vjy) > 0 &&
                 Geometry.triangle_signed_area(vi_1x, vi_1y, vix, viy, vj_1x, vj_1y) <= 0
                p = line_intersect(vi_1x, vi_1y, vix, viy, vjx, vjy, vj_1x, vj_1y)
                vi_p1x, vi_p1y = bayazit_at(i + 1, coords)
                if Geometry.triangle_signed_area(vi_p1x, vi_p1y, vix, viy, p[0], p[1]) < 0
                  d = Geometry.dist_squared(vix, viy, p[0], p[1])
                  if d < lower_dist
                    lower_dist = d
                    lower_int = p
                    lower_index = j
                  end
                end
              end

              vj_p1x, vj_p1y = bayazit_at(j + 1, coords)
              vi_p1x, vi_p1y = bayazit_at(i + 1, coords)
              if Geometry.triangle_signed_area(vi_p1x, vi_p1y, vix, viy, vj_p1x, vj_p1y) > 0 &&
                 Geometry.triangle_signed_area(vi_p1x, vi_p1y, vix, viy, vjx, vjy) <= 0
                p = line_intersect(vi_p1x, vi_p1y, vix, viy, vjx, vjy, vj_p1x, vj_p1y)
                vi_1x, vi_1y = bayazit_at(i - 1, coords)
                if Geometry.triangle_signed_area(vi_1x, vi_1y, vix, viy, p[0], p[1]) > 0
                  d = Geometry.dist_squared(vix, viy, p[0], p[1])
                  if d < upper_dist
                    upper_dist = d
                    upper_index = j
                    upper_int = p
                  end
                end
              end
              j += 1
            end

            if lower_index == (upper_index + 1) % n
              px = (lower_int[0] + upper_int[0]) / 2.0
              py = (lower_int[1] + upper_int[1]) / 2.0

              lower_poly = bayazit_copy(i, upper_index, coords)
              lower_poly << px << py
              upper_poly = bayazit_copy(lower_index, i, coords)
              upper_poly << px << py
            else
              highest_score = 0.0
              best_index = lower_index

              upper_index += n while upper_index < lower_index

              jj = lower_index
              while jj <= upper_index
                if bayazit_can_see?(i, jj, coords)
                  vix,  viy  = bayazit_at(i, coords)
                  vjx,  vjy  = bayazit_at(jj, coords)
                  score = 1.0 / (Geometry.dist_squared(vix, viy, vjx, vjy) + 1)
                  if bayazit_reflex?(jj, coords)
                    vj_1x, vj_1y = bayazit_at(jj - 1, coords)
                    vj_p1x, vj_p1y = bayazit_at(jj + 1, coords)
                    if Geometry.triangle_signed_area(vj_1x, vj_1y, vjx, vjy, vix, viy) <= 0 &&
                       Geometry.triangle_signed_area(vj_p1x, vj_p1y, vjx, vjy, vix, viy) >= 0
                      score += 3
                    else
                      score += 2
                    end
                  else
                    score += 1
                  end
                  if score > highest_score
                    best_index = jj
                    highest_score = score
                  end
                end
                jj += 1
              end
              lower_poly = bayazit_copy(i, best_index, coords)
              upper_poly = bayazit_copy(best_index, i, coords)
            end

            bayazit_decompose(lower_poly, list)
            bayazit_decompose(upper_poly, list)
            return list
          end
          i += 1
        end

        # Polygon is already convex
        list << coords
        list
      end

      # Get vertex [x,y] at index with wrapping. coords: flat array.
      def bayazit_at(i, coords)
        n = coords.length.div(2)
        return [coords[0], coords[1]] if n == 0
        idx = if i < 0
                n - 1 - ((-i - 1) % n)
              else
                i % n
              end
        [coords[idx * 2], coords[idx * 2 + 1]]
      end

      # Copy vertices from index i to j (inclusive) as flat array.
      def bayazit_copy(i, j, coords)
        n = coords.length.div(2)
        j += n while j < i
        result = []
        idx = i
        while idx <= j
          vx, vy = bayazit_at(idx, coords)
          result << vx << vy
          idx += 1
        end
        result
      end

      # Check if vertex at index i can see vertex at index j. coords: flat.
      def bayazit_can_see?(i, j, coords)
        n = coords.length.div(2)
        vix, viy = bayazit_at(i, coords)
        vjx, vjy = bayazit_at(j, coords)

        if bayazit_reflex?(i, coords)
          vi_1x, vi_1y = bayazit_at(i - 1, coords)
          vi_p1x, vi_p1y = bayazit_at(i + 1, coords)
          if Geometry.triangle_signed_area(vix, viy, vi_1x, vi_1y, vjx, vjy) >= 0 &&
             Geometry.triangle_signed_area(vix, viy, vi_p1x, vi_p1y, vjx, vjy) <= 0
            return false
          end
        else
          vi_p1x, vi_p1y = bayazit_at(i + 1, coords)
          vi_1x, vi_1y = bayazit_at(i - 1, coords)
          if Geometry.triangle_signed_area(vix, viy, vi_p1x, vi_p1y, vjx, vjy) <= 0 ||
             Geometry.triangle_signed_area(vix, viy, vi_1x, vi_1y, vjx, vjy) >= 0
            return false
          end
        end

        if bayazit_reflex?(j, coords)
          vj_1x, vj_1y = bayazit_at(j - 1, coords)
          vj_p1x, vj_p1y = bayazit_at(j + 1, coords)
          if Geometry.triangle_signed_area(vjx, vjy, vj_1x, vj_1y, vix, viy) >= 0 &&
             Geometry.triangle_signed_area(vjx, vjy, vj_p1x, vj_p1y, vix, viy) <= 0
            return false
          end
        else
          vj_p1x, vj_p1y = bayazit_at(j + 1, coords)
          vj_1x, vj_1y = bayazit_at(j - 1, coords)
          if Geometry.triangle_signed_area(vjx, vjy, vj_p1x, vj_p1y, vix, viy) <= 0 ||
             Geometry.triangle_signed_area(vjx, vjy, vj_1x, vj_1y, vix, viy) >= 0
            return false
          end
        end

        k = 0
        while k < n
          if !((k + 1) % n == i || k == i || (k + 1) % n == j || k == j)
            vkx, vky = bayazit_at(k, coords)
            vk_p1x, vk_p1y = bayazit_at(k + 1, coords)
            if Geometry.segments_intersect?(vix, viy, vjx, vjy, vkx, vky, vk_p1x, vk_p1y)
              return false
            end
          end
          k += 1
        end
        true
      end

      # Check if vertex at index i is reflex (concave). coords: flat.
      def bayazit_reflex?(i, coords)
        v1x, v1y = bayazit_at(i - 1, coords)
        v2x, v2y = bayazit_at(i, coords)
        v3x, v3y = bayazit_at(i + 1, coords)
        Geometry.triangle_signed_area(v1x, v1y, v2x, v2y, v3x, v3y) < 0
      end

      # Line intersection of (p1→p2) and (p3→p4). Returns [x, y].
      def line_intersect(p1x, p1y, p2x, p2y, p3x, p3y, p4x, p4y)
        a1 = p2y - p1y
        b1 = p1x - p2x
        c1 = a1 * p1x + b1 * p1y

        a2 = p4y - p3y
        b2 = p3x - p4x
        c2 = a2 * p3x + b2 * p3y

        det = a1 * b2 - a2 * b1
        return [0.0, 0.0] if det.abs < 1e-10

        x = (b2 * c1 - b1 * c2) / det
        y = (a1 * c2 - a2 * c1) / det
        [x, y]
      end
    end
  end
end

module Triangles
  # TriVis: Visibility polygon computation by Triangular Expansion
  module TriVis
    class << self

      # Ray-segment intersection returning relative location on ray
      # Returns Float::INFINITY if no intersection
      def seg_intersect_ray(s1x, s1y, s2x, s2y, r1x, r1y, r2x, r2y)
        rdx = (r2x - r1x).to_f
        rdy = (r2y - r1y).to_f

        sdx = (s2x - s1x).to_f
        sdy = (s2y - s1y).to_f

        cross = sdx * rdy - sdy * rdx
        return Float::INFINITY if cross.abs < 1e-10

        t2 = (rdx * (s1y - r1y) + rdy * (r1x - s1x)) / cross

        t1 = if rdx.abs > rdy.abs
               (s1x + sdx * t2 - r1x) / rdx
             else
               (s1y + sdy * t2 - r1y) / rdy
             end

        return Float::INFINITY if t1 < -1e-10 || t2 < -1e-10 || t2 > 1.0 + 1e-10

        t1
      end

      def ray_seg_intersection_point(qx, qy, px, py, s1x, s1y, s2x, s2y)
        t = seg_intersect_ray(s1x, s1y, s2x, s2y, qx, qy, px, py)
        return nil if t == Float::INFINITY
        dx = px - qx
        dy = py - qy
        [qx + t * dx, qy + t * dy]
      end

      def next_edge(e)
        (e % 3 == 2) ? e - 2 : e + 1
      end

      def prev_edge(e)
        (e % 3 == 0) ? e + 2 : e - 1
      end

      def edges_of_tri(t)
        [t * 3, t * 3 + 1, t * 3 + 2]
      end

      def points_of_tri(triangles, t)
        [triangles[t * 3], triangles[t * 3 + 1], triangles[t * 3 + 2]]
      end

      def left_of?(x1, y1, x2, y2, px, py)
        Geometry.triangle_signed_area(x1, y1, x2, y2, px, py) > 0
      end

      def right_of?(x1, y1, x2, y2, px, py)
        Geometry.triangle_signed_area(x1, y1, x2, y2, px, py) < 0
      end

      def orientation(ax, ay, bx, by, cx, cy)
        det = Geometry.triangle_signed_area(ax, ay, bx, by, cx, cy)
        if det > 0
          :ccw
        elsif det < 0
          :cw
        else
          :collinear
        end
      end

      def centroid(coords, triangles, tri)
        pts = points_of_tri(triangles, tri)
        p1x = coords[pts[0] * 2]
        p1y = coords[pts[0] * 2 + 1]
        p2x = coords[pts[1] * 2]
        p2y = coords[pts[1] * 2 + 1]
        p3x = coords[pts[2] * 2]
        p3y = coords[pts[2] * 2 + 1]
        [(p1x + p2x + p3x) / 3.0, (p1y + p2y + p3y) / 3.0]
      end

      def order_angles(qx, qy, p1x, p1y, p2x, p2y)
        seg_left = left_of?(qx, qy, p2x, p2y, p1x, p1y)
        if seg_left
          [p1x, p1y, p2x, p2y]
        else
          [p2x, p2y, p1x, p1y]
        end
      end

      def order_del_angles(coords, qx, qy, p1, p2)
        p1x = coords[p1 * 2]
        p1y = coords[p1 * 2 + 1]
        p2x = coords[p2 * 2]
        p2y = coords[p2 * 2 + 1]
        order_angles(qx, qy, p1x, p1y, p2x, p2y)
      end

      def within_cone?(px, py, slx, sly, srx, sry, rlx, rly, rrx, rry)
        return false if left_of?(px, py, slx, sly, rrx, rry)
        return false if left_of?(px, py, rlx, rly, srx, sry)
        return false if rrx == slx && rry == sly
        return false if srx == rlx && sry == rly
        true
      end

      def restrict_angles(px, py, slx, sly, srx, sry, rlx, rly, rrx, rry)
        nlx = rlx
        nly = rly
        res_left = false
        if right_of?(px, py, rlx, rly, slx, sly)
          nlx = slx
          nly = sly
          res_left = true
        end

        nrx = rrx
        nry = rry
        res_right = false
        if left_of?(px, py, rrx, rry, srx, sry)
          nrx = srx
          nry = sry
          res_right = true
        end

        [nlx, nly, nrx, nry, res_left, res_right]
      end

      # Find triangle containing a point using walk algorithm
      def containing_triangle(coords, triangles, half_edges, qx, qy)
        num_triangles = triangles.length.div(3)
        return -1 if num_triangles == 0

        t = 0
        visited = {}

        iter = 0
        while iter < 100
          return -1 if visited[t]
          visited[t] = true

          e0 = t * 3
          found = true

          ei = 0
          while ei < 3
            e = e0 + ei
            p1 = triangles[e]
            p2 = triangles[next_edge(e)]
            p1x = coords[p1 * 2]
            p1y = coords[p1 * 2 + 1]
            p2x = coords[p2 * 2]
            p2y = coords[p2 * 2 + 1]

            if Geometry.triangle_signed_area(p1x, p1y, p2x, p2y, qx, qy) > 0
              adj = half_edges[e]
              return -1 if adj == -1
              t = adj.div(3)
              found = false
              break
            end
            ei += 1
          end

          return t if found
          iter += 1
        end

        -1
      end

      # Compute visibility polygon from a query point
      def triangular_expansion(del, qx, qy, obstructs = nil, ilx = nil, ily = nil, irx = nil, iry = nil)
        coords = del[:coords]
        triangles = del[:triangles]
        half_edges = del[:half_edges]

        obstructs ||= ->(_) { false }

        tri_start = containing_triangle(coords, triangles, half_edges, qx, qy)
        return [] if tri_start == -1

        cx, cy = centroid(coords, triangles, tri_start)
        start_edges = edges_of_tri(tri_start)
        ret = []
        prestrict = !ilx.nil?

        if prestrict
          ilx, ily, irx, iry = order_angles(qx, qy, ilx, ily, irx, iry)
        end

        i = 0
        while i < 3
          edg = start_edges[i]
          p1 = triangles[edg]
          p2 = triangles[next_edge(edg)]
          slx, sly, srx, sry = order_del_angles(coords, cx, cy, p1, p2)

          return [] if Geometry.dist_squared(slx, sly, qx, qy) == 0.0 || Geometry.dist_squared(srx, sry, qx, qy) == 0.0

          if prestrict
            int_l = seg_intersect_ray(slx, sly, srx, sry, qx, qy, ilx, ily)
            if int_l != Float::INFINITY
              dx = ilx - qx
              dy = ily - qy
              _, _, slx, sly = order_angles(qx, qy, qx + int_l * dx, qy + int_l * dy, slx, sly)
            end
            int_r = seg_intersect_ray(slx, sly, srx, sry, qx, qy, irx, iry)
            if int_r != Float::INFINITY
              dx = irx - qx
              dy = iry - qy
              srx, sry, _, _ = order_angles(qx, qy, srx, sry, qx + int_r * dx, qy + int_r * dy)
            end

            slx, sly, srx, sry = order_angles(qx, qy, slx, sly, srx, sry)
          else
            ilx = slx
            ily = sly
            irx = srx
            iry = sry
          end

          rlx, rly, rrx, rry, _, _ = restrict_angles(qx, qy, slx, sly, srx, sry, ilx, ily, irx, iry)

          if within_cone?(qx, qy, slx, sly, srx, sry, ilx, ily, irx, iry)
            adj = half_edges[edg]
            if adj == -1 || obstructs.call(edg)
              ret << [srx, sry, slx, sly]
            else
              ret.concat(expand(coords, triangles, half_edges, qx, qy, obstructs, adj, rlx, rly, rrx, rry))
            end
          end

          i += 1
        end

        ret
      end

      def expand(coords, triangles, half_edges, qx, qy, obstructs, edg_in, rlx, rly, rrx, rry)
        ret = []
        edg_a = next_edge(edg_in)
        edg_b = prev_edge(edg_in)

        ei = 0
        while ei < 2
          edg = ei == 0 ? edg_a : edg_b
          ei += 1

          p1 = triangles[edg]
          p2 = triangles[next_edge(edg)]
          adj_out = half_edges[edg]

          slx, sly, srx, sry = order_del_angles(coords, qx, qy, p1, p2)

          next unless within_cone?(qx, qy, slx, sly, srx, sry, rlx, rly, rrx, rry)

          nlx, nly, nrx, nry, res_l, res_r = restrict_angles(qx, qy, slx, sly, srx, sry, rlx, rly, rrx, rry)

          next if Geometry.triangle_signed_area(qx, qy, nrx, nry, nlx, nly) <= 0.0

          if adj_out != -1 && !obstructs.call(edg)
            ret.concat(expand(coords, triangles, half_edges, qx, qy, obstructs, adj_out, nlx, nly, nrx, nry))
            next
          end

          unless res_l
            int = seg_intersect_ray(slx, sly, srx, sry, qx, qy, rlx, rly)
            if int != Float::INFINITY
              dx = rlx - qx
              dy = rly - qy
              slx = qx + int * dx
              sly = qy + int * dy
            end
          end

          unless res_r
            int = seg_intersect_ray(slx, sly, srx, sry, qx, qy, rrx, rry)
            if int != Float::INFINITY
              dx = rrx - qx
              dy = rry - qy
              srx = qx + int * dx
              sry = qy + int * dy
            end
          end

          ret << [srx, sry, slx, sly]
        end

        ret
      end

      # Constraint Intersection Handling
      # 1. Find all intersection points between constraint edges
      # 2. Insert intersection points into the point set
      # 3. Split edges at intersection points
      # Call this BEFORE creating the Delaunay triangulation
      def split_intersecting_constraints(coords, edges)
        solver = csolver_init(coords, edges)
        csolver_solve(solver)
      end

      def compute_segment_intersection_params(p1x, p1y, p2x, p2y, p3x, p3y, p4x, p4y)
        d1x = p2x - p1x
        d1y = p2y - p1y
        d2x = p4x - p3x
        d2y = p4y - p3y

        cross = d1x * d2y - d1y * d2x
        return [nil, nil, nil, nil] if cross.abs < 1e-10

        t1 = ((p3x - p1x) * d2y - (p3y - p1y) * d2x) / cross.to_f
        t2 = ((p3x - p1x) * d1y - (p3y - p1y) * d1x) / cross.to_f

        return [nil, nil, nil, nil] unless (0..1).cover?(t1) && (0..1).cover?(t2)

        ix = p1x + t1 * d1x
        iy = p1y + t1 * d1y

        [ix, iy, t1, t2]
      end

      CONSTRAINT_SOLVER_EPSILON = 1e-9 # For floating point comparisons
      CONSTRAINT_SOLVER_SNAP_EPSILON = 1e-6 # For snapping nearby points together

      # Initialize constraint intersection solver
      # @param coords [Array] Flat coordinate array [x1, y1, x2, y2, ...]
      # @param edges [Array] Array of [p1, p2] edge pairs
      # @return [Hash] Solver state
      def csolver_init(coords, edges)
        point_map = {}
        coords_copy = coords.dup

        num_points = coords_copy.length.div(2)
        i = 0
        while i < num_points
          x, y = coords_copy[i * 2], coords_copy[i * 2 + 1]
          key = csolver_point_key(x, y)
          point_map[key] = i
          i += 1
        end

        {
          original_coords: coords.dup,
          original_edges: edges.dup,
          coords: coords_copy,
          point_map: point_map,
          edge_splits: nil
        }
      end

      # Solve constraint intersections
      # @param solver [Hash] Solver state
      # @return [Array] [coords, edges] with intersections resolved
      def csolver_solve(solver)
        csolver_find_and_add_intersections(solver)
        csolver_split_edges_at_interior_points(solver)
        final_edges = csolver_build_split_edges(solver)

        [solver[:coords], final_edges]
      end

      # Generate a point key for lookup
      # @param x [Float] X coordinate
      # @param y [Float] Y coordinate
      # @return [String] Point key
      def csolver_point_key(x, y)
        "#{x.round(8)},#{y.round(8)}"
      end

      # Get or create a point, snapping to nearby existing points
      # @param solver [Hash] Solver state
      # @param x [Float] X coordinate
      # @param y [Float] Y coordinate
      # @return [Integer] Point index
      def csolver_get_or_create_point(solver, x, y)
        coords = solver[:coords]
        point_map = solver[:point_map]
        key = csolver_point_key(x, y)

        return point_map[key] if point_map[key]

        point_map.each do |_k, idx|
          px, py = coords[idx * 2], coords[idx * 2 + 1]
          if (px - x).abs < CONSTRAINT_SOLVER_SNAP_EPSILON &&
             (py - y).abs < CONSTRAINT_SOLVER_SNAP_EPSILON
            return idx
          end
        end

        new_idx = coords.length.div(2)
        coords << x << y
        point_map[key] = new_idx
        new_idx
      end

      # Get point coordinates
      # @param solver [Hash] Solver state
      # @param idx [Integer] Point index
      # @return [Array] [x, y]
      def csolver_point_coords(solver, idx)
        coords = solver[:coords]
        [coords[idx * 2], coords[idx * 2 + 1]]
      end

      # Check if two segments share an endpoint
      # @param e1 [Array] [p1, p2]
      # @param e2 [Array] [p3, p4]
      # @return [Boolean]
      def csolver_segments_share_endpoint?(e1, e2)
        p1, p2 = e1
        p3, p4 = e2
        p1 == p3 || p1 == p4 || p2 == p3 || p2 == p4
      end

      # Compute intersection between two segments
      # @param solver [Hash] Solver state
      # @param e1 [Array] [p1, p2]
      # @param e2 [Array] [p3, p4]
      # @return [Array] [ix, iy, t1, t2] or [nil, nil, nil, nil]
      def csolver_segment_intersection(solver, e1, e2)
        p1, p2 = e1
        p3, p4 = e2

        p1x, p1y = csolver_point_coords(solver, p1)
        p2x, p2y = csolver_point_coords(solver, p2)
        p3x, p3y = csolver_point_coords(solver, p3)
        p4x, p4y = csolver_point_coords(solver, p4)

        compute_segment_intersection_params(p1x, p1y, p2x, p2y, p3x, p3y, p4x, p4y)
      end

      # Check if parameter t represents an interior intersection
      # @param t [Float] Parameter value
      # @param epsilon [Float] Tolerance
      # @return [Boolean]
      def csolver_interior_intersection?(t, epsilon = 0.001)
        t > epsilon && t < (1.0 - epsilon)
      end

      # Find all intersections and add them as new points
      # @param solver [Hash] Solver state
      def csolver_find_and_add_intersections(solver)
        original_edges = solver[:original_edges]
        max_iterations = 100
        iteration = 0

        found_intersection = true
        while found_intersection && iteration < max_iterations
          iteration += 1
          found_intersection = false

          i = 0
          while i < original_edges.length
            e1 = original_edges[i]
            j = i + 1
            while j < original_edges.length
              e2 = original_edges[j]
              unless csolver_segments_share_endpoint?(e1, e2)
                ix, iy, t1, t2 = csolver_segment_intersection(solver, e1, e2)

                if ix && csolver_interior_intersection?(t1) &&
                   csolver_interior_intersection?(t2)
                  existing_count = solver[:coords].length.div(2)
                  new_idx = csolver_get_or_create_point(solver, ix, iy)

                  if new_idx >= existing_count
                    found_intersection = true
                  end
                end
              end
              j += 1
            end
            i += 1
          end
        end
      end

      # Check if the point is on the segment interior
      # @param px [Float] Point X
      # @param py [Float] Point Y
      # @param p1x [Float] Segment start X
      # @param p1y [Float] Segment start Y
      # @param p2x [Float] Segment end X
      # @param p2y [Float] Segment end Y
      # @return [Boolean]
      def csolver_point_on_segment_interior?(px, py, p1x, p1y, p2x, p2y)
        cross = (p2x - p1x) * (py - p1y) - (p2y - p1y) * (px - p1x)
        return false if cross.abs > CONSTRAINT_SOLVER_EPSILON

        if (p2x - p1x).abs > CONSTRAINT_SOLVER_EPSILON
          t = (px - p1x) / (p2x - p1x)
        elsif (p2y - p1y).abs > CONSTRAINT_SOLVER_EPSILON
          t = (py - p1y) / (p2y - p1y)
        else
          return false
        end

        t > 0.001 && t < 0.999
      end

      # Get parameter value for point on a segment
      # @param px [Float] Point X
      # @param py [Float] Point Y
      # @param p1x [Float] Segment start X
      # @param p1y [Float] Segment start Y
      # @param p2x [Float] Segment end X
      # @param p2y [Float] Segment end Y
      # @return [Float] Parameter t
      def csolver_point_parameter_on_segment(px, py, p1x, p1y, p2x, p2y)
        dx = p2x - p1x
        dy = p2y - p1y

        if dx.abs > dy.abs
          (px - p1x) / dx
        elsif dy.abs > CONSTRAINT_SOLVER_EPSILON
          (py - p1y) / dy
        else
          0.0
        end
      end

      # Find and record all points that lie on edge interiors
      # @param solver [Hash] Solver state
      def csolver_split_edges_at_interior_points(solver)
        edge_splits = {}
        original_edges = solver[:original_edges]
        coords = solver[:coords]

        edge_idx = 0
        while edge_idx < original_edges.length
          edge = original_edges[edge_idx]
          p1, p2 = edge
          p1x, p1y = csolver_point_coords(solver, p1)
          p2x, p2y = csolver_point_coords(solver, p2)

          splits = []

          num_points = coords.length.div(2)
          pt_idx = 0
          while pt_idx < num_points
            if pt_idx != p1 && pt_idx != p2
              px, py = csolver_point_coords(solver, pt_idx)

              if csolver_point_on_segment_interior?(px, py, p1x, p1y, p2x, p2y)
                t = csolver_point_parameter_on_segment(px, py, p1x, p1y, p2x, p2y)
                splits << { idx: pt_idx, t: t }
              end
            end
            pt_idx += 1
          end

          edge_splits[edge_idx] = splits.sort_by { |s| s[:t] } unless splits.empty?
          edge_idx += 1
        end

        solver[:edge_splits] = edge_splits
      end

      # Build the final edge list with splits applied
      # @param solver [Hash] Solver state
      # @return [Array] Array of [p1, p2] edges
      def csolver_build_split_edges(solver)
        result = []
        original_edges = solver[:original_edges]
        edge_splits = solver[:edge_splits]

        edge_idx = 0
        while edge_idx < original_edges.length
          edge = original_edges[edge_idx]
          p1, p2 = edge
          splits = edge_splits[edge_idx]

          if splits.nil? || splits.empty?
            result << edge
          else
            current = p1
            si = 0
            while si < splits.length
              result << [current, splits[si][:idx]]
              current = splits[si][:idx]
              si += 1
            end
            result << [current, p2]
          end
          edge_idx += 1
        end

        result
      end

      # Initialize constrained Delaunay triangulation from a Delaunay result.
      # del: {coords:, triangles:, half_edges:}. edges: optional [[p1,p2], ...].
      # Returns constraint state hash.
      def cdt_init(del, edges = nil)
        num_points = del[:coords].length.div(2)
        num_edges = del[:triangles].length

        vert_map = Array.new(num_points)
        flips = {}
        consd = {}

        con = {
          del: del,
          vert_map: vert_map,
          flips: flips,
          consd: consd
        }

        e = 0
        while e < num_edges
          v = del[:triangles][e]
          cdt_update_vert(con, e) if vert_map[v].nil?
          e += 1
        end

        # Pre-delaunify: converge to fully Delaunay before constraining.
        # Converge to fully Delaunay before constraining edges.
        cdt_delaunify(con, deep: true)

        cdt_constrain_edges(con, edges) if edges
        con
      end

      # Restore Delaunay condition for all non-constrained edges.
      # deep: if true, repeat until fully converged.
      # Uses inline flips for speed, then rebuilds vert_map once at the end.
      def cdt_delaunify(con, deep: false)
        del = con[:del]
        triangles = del[:triangles]
        half_edges = del[:half_edges]
        consd = con[:consd]
        num_edges = triangles.length
        max_passes = deep ? num_edges : 1

        pass = 0
        keep_going = true
        while keep_going
          flipped = 0
          edg = 0
          while edg < num_edges
            unless consd[edg]
              adj = half_edges[edg]
              unless adj == -1 || cdt_delaunay_fast?(con, edg)
                bot = prev_edge(edg)
                top = prev_edge(adj)
                adj_bot = half_edges[bot]
                adj_top = half_edges[top]

                triangles[edg] = triangles[top]
                half_edges[edg] = adj_top
                half_edges[adj_top] = edg if adj_top != -1
                half_edges[bot] = top

                triangles[adj] = triangles[bot]
                half_edges[adj] = adj_bot
                half_edges[adj_bot] = adj if adj_bot != -1
                half_edges[top] = bot

                flipped += 1
              end
            end
            edg += 1
          end
          pass += 1
          keep_going = deep && flipped > 0 && pass < max_passes
        end

        # Rebuild vert_map from scratch — faster than per-flip updates and
        # guarantees correctness after all flips are done.
        vert_map = con[:vert_map]
        i = 0
        while i < vert_map.length
          vert_map[i] = nil
          i += 1
        end
        e = 0
        while e < num_edges
          v = triangles[e]
          cdt_update_vert(con, e) if vert_map[v].nil?
          e += 1
        end

        con[:flips].clear
      end

      # Constrain a single edge between two point indices. Returns edge id.
      def cdt_constrain_edge(con, seg_p1, seg_p2)
        del = con[:del]
        vert_map = con[:vert_map]
        consd = con[:consd]
        triangles = del[:triangles]
        half_edges = del[:half_edges]
        start = vert_map[seg_p1]

        # Walk the edges touching seg_p1
        edg = start
        while edg != -1
          # edg points toward seg_p1, so its start-point is opposite it
          p4 = triangles[edg]
          nxt = next_edge(edg)

          # Already constrained, but in reverse order
          if p4 == seg_p2
            return cdt_protect(con, edg)
          end

          # The edge opposite seg_p1
          opp = prev_edge(edg)
          p3 = triangles[opp]

          # Already constrained
          if p3 == seg_p2
            cdt_protect(con, nxt)
            return nxt
          end

          # Edge opposite seg_p1 intersects constraint
          if cdt_intersect_segments?(con, seg_p1, seg_p2, p3, p4)
            edg = opp
            break
          end

          adj = half_edges[nxt]
          edg = adj
          break if edg == start
        end

        con_edge = edg
        rescan = -1

        # Walk through triangulation, flipping intersecting edges
        while edg != -1
          adj = half_edges[edg]
          bot = prev_edge(edg)
          top = prev_edge(adj)
          rgt = next_edge(adj)

          raise "Constraining edge exited the hull" if adj == -1
          raise "Edge intersects already constrained edge" if consd[edg]

          if cdt_collinear?(con, seg_p1, seg_p2, triangles[edg]) ||
             cdt_collinear?(con, seg_p1, seg_p2, triangles[adj])
            raise "Constraining edge intersects point"
          end

          # Convexity test: the quadrilateral is convex if its diagonals intersect
          convex = cdt_intersect_segments?(
            con,
            triangles[edg], triangles[adj],
            triangles[bot], triangles[top]
          )

          unless convex
            # Non-convex: can't flip, continue and rescan later
            rescan = edg if rescan == -1

            if triangles[top] == seg_p2
              raise "Infinite loop: non-convex quadrilateral" if edg == rescan
              edg = rescan
              rescan = -1
              next
            end

            # Look for the next intersection
            if cdt_intersect_segments?(con, seg_p1, seg_p2, triangles[top], triangles[adj])
              edg = top
            elsif cdt_intersect_segments?(con, seg_p1, seg_p2, triangles[rgt], triangles[top])
              edg = rgt
            elsif rescan == edg
              raise "Infinite loop: no further intersect after non-convex"
            end

            next
          end

          cdt_flip_diagonal(con, edg)

          # New edge might still intersect
          if cdt_intersect_segments?(con, seg_p1, seg_p2, triangles[bot], triangles[top])
            rescan = bot if rescan == -1
            raise "Infinite loop: flipped diagonal still intersects" if rescan == bot
          end

          # Reached the other endpoint?
          if triangles[top] == seg_p2
            con_edge = top
            edg = rescan
            rescan = -1
          elsif cdt_intersect_segments?(con, seg_p1, seg_p2, triangles[rgt], triangles[top])
            edg = rgt
          end
          # Note: if neither condition is true, edg is unchanged and we continue
          # with the same edge (which now has different vertex configuration after flip)
        end

        cdt_protect(con, con_edge)

        # Restore Delaunay condition for flipped edges.
        # The snapshot flips each pass, deletes both edg+adj,
        # repeat while any edge was flipped. New entries from cdt_flip_diagonal
        # are picked up in the next pass.
        flips = con[:flips]
        flipped = 1
        while flipped > 0
          flipped = 0
          keys = flips.keys
          ki = 0
          while ki < keys.length
            e = keys[ki]
            ki += 1
            next unless flips.delete(e)  # skip if already removed as adj
            adj_e = half_edges[e]
            next if adj_e == -1
            flips.delete(adj_e)
            unless cdt_delaunay_fast?(con, e)
              cdt_flip_diagonal(con, e)
              flipped += 1
            end
          end
        end

        cdt_find_edge(con, seg_p1, seg_p2)
      end

      # Constrain all edges
      # @param con [Hash] Constraint state
      # @param edges [Array] Array of [p1, p2] pairs
      # @return [Hash] Constraint state
      def cdt_constrain_edges(con, edges)
        i = 0
        while i < edges.length
          e = edges[i]
          cdt_constrain_edge(con, e[0], e[1])
          i += 1
        end
        con
      end

      # Check if edge is constrained
      # @param con [Hash] Constraint state
      # @param edg [Integer] Edge id
      # @return [Boolean]
      def cdt_edge_constrained?(con, edg)
        con[:consd][edg] || false
      end

      # Find edge from p1 to p2
      # @param con [Hash] Constraint state
      # @param p1 [Integer] First point index
      # @param p2 [Integer] Second point index
      # @return [Integer] Edge id, negative if on hull pointing other way, or nil
      def cdt_find_edge(con, p1, p2)
        vert_map = con[:vert_map]
        del = con[:del]
        start1 = vert_map[p2]
        triangles = del[:triangles]
        half_edges = del[:half_edges]

        edg = start1
        prv = -1
        while edg != -1
          return edg if triangles[edg] == p1
          prv = next_edge(edg)
          edg = half_edges[prv]
          break if edg == start1
        end

        # Check hull edge
        return -prv if triangles[next_edge(prv)] == p1
        nil
      end

      # Mark edge as constrained
      # @param con [Hash] Constraint state
      # @param edg [Integer] Edge id
      # @return [Integer] Adjacent edge id or negative edge id
      def cdt_protect(con, edg)
        del = con[:del]
        flips = con[:flips]
        consd = con[:consd]
        adj = del[:half_edges][edg]
        flips.delete(edg)
        consd[edg] = true

        if adj != -1
          flips.delete(adj)
          consd[adj] = true
          return adj
        end

        -edg
      end

      # Mark the edge as flipped (unless constrained)
      # @param con [Hash] Constraint state
      # @param edg [Integer] Edge id
      # @return [Boolean] True if marked, false if constrained
      def cdt_mark_flip(con, edg)
        consd = con[:consd]
        return false if consd[edg]
        flips = con[:flips]
        del = con[:del]
        adj = del[:half_edges][edg]
        if adj != -1
          flips[edg] = true
          flips[adj] = true
        end
        true
      end

      # Flip the diagonal of two adjacent triangles
      # @param con [Hash] Constraint state
      # @param edg [Integer] Edge id
      # @return [Integer] New bottom edge id
      def cdt_flip_diagonal(con, edg)
        del = con[:del]
        flips = con[:flips]
        consd = con[:consd]
        triangles = del[:triangles]
        half_edges = del[:half_edges]

        adj = half_edges[edg]
        bot = prev_edge(edg)
        lft = next_edge(edg)
        top = prev_edge(adj)
        rgt = next_edge(adj)
        adj_bot = half_edges[bot]
        adj_top = half_edges[top]

        raise "Trying to flip constrained edge" if consd[edg]

        # Move edg to top: copy top's vertex and connectivity to edg
        triangles[edg] = triangles[top]
        half_edges[edg] = adj_top
        # Propagate top's flip/constrained status to edg
        if flips[top]
          flips[edg] = true
        else
          flips.delete(edg)
          # Only propagate constrained status when not flipped
          if consd[top]
            consd[edg] = true
          else
            consd.delete(edg)
          end
        end
        half_edges[adj_top] = edg if adj_top != -1
        half_edges[bot] = top

        # Move adj to bot: copy bot's vertex and connectivity to adj
        triangles[adj] = triangles[bot]
        half_edges[adj] = adj_bot
        # Propagate bot's flip/constrained status to adj
        if flips[bot]
          flips[adj] = true
        else
          flips.delete(adj)
          # Only propagate constrained status when not flipped
          if consd[bot]
            consd[adj] = true
          else
            consd.delete(adj)
          end
        end
        half_edges[adj_bot] = adj if adj_bot != -1
        half_edges[top] = bot

        cdt_mark_flip(con, edg)
        cdt_mark_flip(con, lft)
        cdt_mark_flip(con, adj)
        cdt_mark_flip(con, rgt)

        # Mark bot and top as flipped unconditionally (they are the new diagonal)
        flips[bot] = true
        consd.delete(bot)
        flips[top] = true
        consd.delete(top)

        cdt_update_vert(con, edg)
        cdt_update_vert(con, lft)
        cdt_update_vert(con, adj)
        cdt_update_vert(con, rgt)

        bot
      end

      # Check if the edge satisfies Delaunay condition
      # Constrained edges are always considered Delaunay (they shouldn't be flipped)
      # @param con [Hash] Constraint state
      # @param edg [Integer] Edge id
      # @return [Boolean]
      def cdt_delaunay?(con, edg)
        consd = con[:consd]
        return true if consd[edg]

        del = con[:del]
        triangles = del[:triangles]
        half_edges = del[:half_edges]
        adj = half_edges[edg]
        return true if adj == -1

        p1 = triangles[prev_edge(edg)]
        p2 = triangles[edg]
        p3 = triangles[next_edge(edg)]
        px = triangles[prev_edge(adj)]

        !cdt_in_circle?(con, p1, p2, p3, px)
      end

      # Delaunay check using robust incircle predicate.
      # Returns true when edge satisfies Delaunay condition (opposite vertex
      # is outside or on circumcircle). Constrained and hull edges always pass.
      def cdt_delaunay_fast?(con, edg)
        consd = con[:consd]
        return true if consd[edg]

        del = con[:del]
        triangles = del[:triangles]
        half_edges = del[:half_edges]
        adj = half_edges[edg]
        return true if adj == -1

        p1 = triangles[prev_edge(edg)]
        p2 = triangles[edg]
        p3 = triangles[next_edge(edg)]
        px = triangles[prev_edge(adj)]

        coords = del[:coords]
        Geometry.robust_incircle(
          coords[p1 * 2], coords[p1 * 2 + 1],
          coords[p2 * 2], coords[p2 * 2 + 1],
          coords[p3 * 2], coords[p3 * 2 + 1],
          coords[px * 2], coords[px * 2 + 1]
        ) >= 0
      end

      # Update vertex -> edge map
      # @param con [Hash] Constraint state
      # @param start [Integer] Starting edge
      # @return [Integer] Updated incoming edge
      def cdt_update_vert(con, start)
        del = con[:del]
        vert_map = con[:vert_map]
        triangles = del[:triangles]
        half_edges = del[:half_edges]
        v = triangles[start]

        # Walk counter-clockwise to find right-most incoming edge
        inc = prev_edge(start)
        adj = half_edges[inc]
        while adj != -1 && adj != start
          inc = prev_edge(adj)
          adj = half_edges[inc]
        end

        vert_map[v] = inc
        inc
      end

      # Check if segments [p1,p2] and [p3,p4] intersect
      # @param con [Hash] Constraint state
      # @param p1 [Integer] First point of first segment
      # @param p2 [Integer] Second point of first segment
      # @param p3 [Integer] First point of second segment
      # @param p4 [Integer] Second point of second segment
      # @return [Boolean]
      def cdt_intersect_segments?(con, p1, p2, p3, p4)
        return false if p1 == p3 || p1 == p4 || p2 == p3 || p2 == p4

        coords = con[:del][:coords]
        cdt_intersect_segments_coords(
          coords[p1 * 2], coords[p1 * 2 + 1],
          coords[p2 * 2], coords[p2 * 2 + 1],
          coords[p3 * 2], coords[p3 * 2 + 1],
          coords[p4 * 2], coords[p4 * 2 + 1]
        )
      end

      # Check segment intersection with coordinates.
      # Port of Constrainautor.ts intersectSegments using robust orient2d.
      # @return [Boolean]
      def cdt_intersect_segments_coords(p1x, p1y, p2x, p2y, p3x, p3y, p4x, p4y)
        x0 = Geometry.robust_orient2d(p1x, p1y, p3x, p3y, p4x, p4y)
        y0 = Geometry.robust_orient2d(p2x, p2y, p3x, p3y, p4x, p4y)
        return false if (x0 > 0 && y0 > 0) || (x0 < 0 && y0 < 0)

        x1 = Geometry.robust_orient2d(p3x, p3y, p1x, p1y, p2x, p2y)
        y1 = Geometry.robust_orient2d(p4x, p4y, p1x, p1y, p2x, p2y)
        return false if (x1 > 0 && y1 > 0) || (x1 < 0 && y1 < 0)

        # Degenerate collinear case
        if x0 == 0 && y0 == 0 && x1 == 0 && y1 == 0
          return !([p3x, p4x].max < [p1x, p2x].min ||
            [p1x, p2x].max < [p3x, p4x].min ||
            [p3y, p4y].max < [p1y, p2y].min ||
            [p1y, p2y].max < [p3y, p4y].min)
        end

        true
      end

      # Check if point px is in circumcircle of CCW (y-down) triangle (p1,p2,p3).
      def cdt_in_circle?(con, p1, p2, p3, px)
        coords = con[:del][:coords]
        Geometry.robust_incircle(
          coords[p1 * 2], coords[p1 * 2 + 1],
          coords[p2 * 2], coords[p2 * 2 + 1],
          coords[p3 * 2], coords[p3 * 2 + 1],
          coords[px * 2], coords[px * 2 + 1]
        ) < 0
      end

      # Check if points p1, p2, p are collinear.
      # @param con [Hash] Constraint state
      # @return [Boolean]
      def cdt_collinear?(con, p1, p2, p)
        coords = con[:del][:coords]
        Geometry.robust_orient2d(
          coords[p1 * 2], coords[p1 * 2 + 1],
          coords[p2 * 2], coords[p2 * 2 + 1],
          coords[p * 2], coords[p * 2 + 1]
        ) == 0.0
      end
    end
  end
end

