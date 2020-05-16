package opensimplex

import kotlin.experimental.xor

/**
 * Converted to kotlin from https://github.com/KdotJPG/OpenSimplex2
 *
 * K.jpg's OpenSimplex 2, smooth variant ("SuperSimplex")
 *
 * - 2D is standard simplex, modified to support larger kernels.
 * Implemented using a lookup table.
 * - 3D is "Re-oriented 8-point BCC noise" which constructs an
 * isomorphic BCC lattice in a much different way than usual.
 * - 4D uses a naïve pregenerated lookup table, and averages out
 * to the expected performance.
 *
 * Multiple versions of each function are provided. See the
 * documentation above each, for more info.
 */
class OpenSimplex2S(seed: Long) {
    private val perm: ShortArray
    private val permGrad2: Array<Grad2?>
    private val permGrad3: Array<Grad3?>
    private val permGrad4: Array<Grad4?>
    /*
	 * Noise Evaluators
	 */
    /**
     * 2D SuperSimplex noise, standard lattice orientation.
     */
    fun noise2(x: Double, y: Double): Double {

        // Get points for A2* lattice
        val s = 0.366025403784439 * (x + y)
        val xs = x + s
        val ys = y + s
        return noise2_Base(xs, ys)
    }

    /**
     * 2D SuperSimplex noise, with Y pointing down the main diagonal.
     * Might be better for a 2D sandbox style game, where Y is vertical.
     * Probably slightly less optimal for heightmaps or continent maps.
     */
    fun noise2_XBeforeY(x: Double, y: Double): Double {

        // Skew transform and rotation baked into one.
        val xx = x * 0.7071067811865476
        val yy = y * 1.224744871380249
        return noise2_Base(yy + xx, yy - xx)
    }

    /**
     * 2D SuperSimplex noise base.
     * Lookup table implementation inspired by DigitalShadow.
     */
    private fun noise2_Base(xs: Double, ys: Double): Double {
        var value = 0.0

        // Get base points and offsets
        val xsb = fastFloor(xs)
        val ysb = fastFloor(ys)
        val xsi = xs - xsb
        val ysi = ys - ysb

        // Index to point list
        val a = (xsi + ysi).toInt()
        val index = a shl 2 or (
                (xsi - ysi / 2 + 1 - a / 2.0).toInt() shl 3) or (
                (ysi - xsi / 2 + 1 - a / 2.0).toInt() shl 4)
        val ssi = (xsi + ysi) * -0.211324865405187
        val xi = xsi + ssi
        val yi = ysi + ssi

        // Point contributions
        for (i in 0..3) {
            val c = LOOKUP_2D[index + i]
            val dx = xi + c!!.dx
            val dy = yi + c.dy
            var attn = 2.0 / 3.0 - dx * dx - dy * dy
            if (attn <= 0) continue
            val pxm = xsb + c.xsv and PMASK
            val pym = ysb + c.ysv and PMASK
            val grad = permGrad2[(perm[pxm] xor pym.toShort()).toInt()]
            val extrapolation = grad!!.dx * dx + grad.dy * dy
            attn *= attn
            value += attn * attn * extrapolation
        }
        return value
    }

    /**
     * 3D Re-oriented 8-point BCC noise, classic orientation
     * Proper substitute for what 3D SuperSimplex would be,
     * in light of Forbidden Formulae.
     * Use noise3_XYBeforeZ or noise3_XZBeforeY instead, wherever appropriate.
     */
    fun noise3_Classic(x: Double, y: Double, z: Double): Double {

        // Re-orient the cubic lattices via rotation, to produce the expected look on cardinal planar slices.
        // If texturing objects that don't tend to have cardinal plane faces, you could even remove this.
        // Orthonormal rotation. Not a skew transform.
        val r = 2.0 / 3.0 * (x + y + z)
        val xr = r - x
        val yr = r - y
        val zr = r - z

        // Evaluate both lattices to form a BCC lattice.
        return noise3_BCC(xr, yr, zr)
    }

    /**
     * 3D Re-oriented 8-point BCC noise, with better visual isotropy in (X, Y).
     * Recommended for 3D terrain and time-varied animations.
     * The Z coordinate should always be the "different" coordinate in your use case.
     * If Y is vertical in world coordinates, call noise3_XYBeforeZ(x, z, Y) or use noise3_XZBeforeY.
     * If Z is vertical in world coordinates, call noise3_XYBeforeZ(x, y, Z).
     * For a time varied animation, call noise3_XYBeforeZ(x, y, T).
     */
    fun noise3_XYBeforeZ(x: Double, y: Double, z: Double): Double {

        // Re-orient the cubic lattices without skewing, to make X and Y triangular like 2D.
        // Orthonormal rotation. Not a skew transform.
        val xy = x + y
        val s2 = xy * -0.211324865405187
        val zz = z * 0.577350269189626
        val xr = x + s2 - zz
        val yr = y + s2 - zz
        val zr = xy * 0.577350269189626 + zz

        // Evaluate both lattices to form a BCC lattice.
        return noise3_BCC(xr, yr, zr)
    }

    /**
     * 3D Re-oriented 8-point BCC noise, with better visual isotropy in (X, Z).
     * Recommended for 3D terrain and time-varied animations.
     * The Y coordinate should always be the "different" coordinate in your use case.
     * If Y is vertical in world coordinates, call noise3_XZBeforeY(x, Y, z).
     * If Z is vertical in world coordinates, call noise3_XZBeforeY(x, Z, y) or use noise3_XYBeforeZ.
     * For a time varied animation, call noise3_XZBeforeY(x, T, y) or use noise3_XYBeforeZ.
     */
    fun noise3_XZBeforeY(x: Double, y: Double, z: Double): Double {

        // Re-orient the cubic lattices without skewing, to make X and Z triangular like 2D.
        // Orthonormal rotation. Not a skew transform.
        val xz = x + z
        val s2 = xz * -0.211324865405187
        val yy = y * 0.577350269189626
        val xr = x + s2 - yy
        val zr = z + s2 - yy
        val yr = xz * 0.577350269189626 + yy

        // Evaluate both lattices to form a BCC lattice.
        return noise3_BCC(xr, yr, zr)
    }

    /**
     * Generate overlapping cubic lattices for 3D Re-oriented BCC noise.
     * Lookup table implementation inspired by DigitalShadow.
     * It was actually faster to narrow down the points in the loop itself,
     * than to build up the index with enough info to isolate 8 points.
     */
    private fun noise3_BCC(xr: Double, yr: Double, zr: Double): Double {

        // Get base and offsets inside cube of first lattice.
        val xrb = fastFloor(xr)
        val yrb = fastFloor(yr)
        val zrb = fastFloor(zr)
        val xri = xr - xrb
        val yri = yr - yrb
        val zri = zr - zrb

        // Identify which octant of the cube we're in. This determines which cell
        // in the other cubic lattice we're in, and also narrows down one point on each.
        val xht = (xri + 0.5).toInt()
        val yht = (yri + 0.5).toInt()
        val zht = (zri + 0.5).toInt()
        val index = xht shl 0 or (yht shl 1) or (zht shl 2)

        // Point contributions
        var value = 0.0
        var c = LOOKUP_3D[index]
        while (c != null) {
            val dxr = xri + c.dxr
            val dyr = yri + c.dyr
            val dzr = zri + c.dzr
            var attn = 0.75 - dxr * dxr - dyr * dyr - dzr * dzr
            if (attn < 0) {
                c = c.nextOnFailure
            } else {
                val pxm = xrb + c.xrv and PMASK
                val pym = yrb + c.yrv and PMASK
                val pzm = zrb + c.zrv and PMASK
                val grad = permGrad3[(perm[(perm[pxm] xor pym.toShort()).toInt()] xor pzm.toShort()).toInt()]
                val extrapolation = grad!!.dx * dxr + grad.dy * dyr + grad.dz * dzr
                attn *= attn
                value += attn * attn * extrapolation
                c = c.nextOnSuccess
            }
        }
        return value
    }

    /**
     * 4D SuperSimplex noise, classic lattice orientation.
     */
    fun noise4_Classic(x: Double, y: Double, z: Double, w: Double): Double {

        // Get points for A4 lattice
        val s = 0.309016994374947 * (x + y + z + w)
        val xs = x + s
        val ys = y + s
        val zs = z + s
        val ws = w + s
        return noise4_Base(xs, ys, zs, ws)
    }

    /**
     * 4D SuperSimplex noise, with XY and ZW forming orthogonal triangular-based planes.
     * Recommended for 3D terrain, where X and Y (or Z and W) are horizontal.
     * Recommended for noise(x, y, sin(time), cos(time)) trick.
     */
    fun noise4_XYBeforeZW(x: Double, y: Double, z: Double, w: Double): Double {
        val s2 = (x + y) * -0.28522513987434876941 + (z + w) * 0.83897065470611435718
        val t2 = (z + w) * 0.21939749883706435719 + (x + y) * -0.48214856493302476942
        val xs = x + s2
        val ys = y + s2
        val zs = z + t2
        val ws = w + t2
        return noise4_Base(xs, ys, zs, ws)
    }

    /**
     * 4D SuperSimplex noise, with XZ and YW forming orthogonal triangular-based planes.
     * Recommended for 3D terrain, where X and Z (or Y and W) are horizontal.
     */
    fun noise4_XZBeforeYW(x: Double, y: Double, z: Double, w: Double): Double {
        val s2 = (x + z) * -0.28522513987434876941 + (y + w) * 0.83897065470611435718
        val t2 = (y + w) * 0.21939749883706435719 + (x + z) * -0.48214856493302476942
        val xs = x + s2
        val ys = y + t2
        val zs = z + s2
        val ws = w + t2
        return noise4_Base(xs, ys, zs, ws)
    }

    /**
     * 4D SuperSimplex noise, with XYZ oriented like noise3_Classic,
     * and W for an extra degree of freedom.
     * Recommended for time-varied animations which texture a 3D object (W=time)
     */
    fun noise4_XYZBeforeW(x: Double, y: Double, z: Double, w: Double): Double {
        val xyz = x + y + z
        val ww = w * 1.118033988749894
        val s2 = xyz / -6.0 + ww
        val xs = x + s2
        val ys = y + s2
        val zs = z + s2
        val ws = -0.5 * xyz + ww
        return noise4_Base(xs, ys, zs, ws)
    }

    /**
     * 4D SuperSimplex noise base.
     * Using ultra-simple 4x4x4x4 lookup partitioning.
     */
    private fun noise4_Base(xs: Double, ys: Double, zs: Double, ws: Double): Double {
        var value = 0.0

        // Get base points and offsets
        val xsb = fastFloor(xs)
        val ysb = fastFloor(ys)
        val zsb = fastFloor(zs)
        val wsb = fastFloor(ws)
        val xsi = xs - xsb
        val ysi = ys - ysb
        val zsi = zs - zsb
        val wsi = ws - wsb

        // Unskewed offsets
        val ssi = (xsi + ysi + zsi + wsi) * -0.138196601125011
        val xi = xsi + ssi
        val yi = ysi + ssi
        val zi = zsi + ssi
        val wi = wsi + ssi
        val index = (fastFloor(xs * 4) and 3 shl 0
                or (fastFloor(ys * 4) and 3 shl 2)
                or (fastFloor(zs * 4) and 3 shl 4)
                or (fastFloor(ws * 4) and 3 shl 6))

        // Point contributions
        for (c in LOOKUP_4D[index]!!) {
            val dx = xi + c!!.dx
            val dy = yi + c.dy
            val dz = zi + c.dz
            val dw = wi + c.dw
            var attn = 0.8 - dx * dx - dy * dy - dz * dz - dw * dw
            if (attn > 0) {
                attn *= attn
                val pxm = xsb + c.xsv and PMASK
                val pym = ysb + c.ysv and PMASK
                val pzm = zsb + c.zsv and PMASK
                val pwm = wsb + c.wsv and PMASK
                val grad =
                    permGrad4[(perm[(perm[(perm[pxm] xor pym.toShort()).toInt()] xor pzm.toShort()).toInt()] xor pwm.toShort()).toInt()]
                val extrapolation = grad!!.dx * dx + grad.dy * dy + grad.dz * dz + grad.dw * dw
                value += attn * attn * extrapolation
            }
        }
        return value
    }

    companion object {
        private const val PSIZE = 2048
        private const val PMASK = 2047

        /*
	 * Utility
	 */
        private fun fastFloor(x: Double): Int {
            val xi = x.toInt()
            return if (x < xi) xi - 1 else xi
        }

        /*
	 * Definitions
	 */
        private val LOOKUP_2D: Array<LatticePoint2D?>
        private val LOOKUP_3D: Array<LatticePoint3D?>
        private val LOOKUP_4D: Array<Array<LatticePoint4D?>?>
        const val N2 = 0.05481866495625118
        const val N3 = 0.2781926117527186
        const val N4 = 0.11127401889945551
        private val GRADIENTS_2D: Array<Grad2?>
        private val GRADIENTS_3D: Array<Grad3?>
        private val GRADIENTS_4D: Array<Grad4?>

        init {
            LOOKUP_2D = arrayOfNulls(8 * 4)
            LOOKUP_3D = arrayOfNulls(8)
            LOOKUP_4D = arrayOfNulls(256)
            for (i in 0..7) {
                var i1: Int
                var j1: Int
                var i2: Int
                var j2: Int
                if (i and 1 == 0) {
                    if (i and 2 == 0) {
                        i1 = -1
                        j1 = 0
                    } else {
                        i1 = 1
                        j1 = 0
                    }
                    if (i and 4 == 0) {
                        i2 = 0
                        j2 = -1
                    } else {
                        i2 = 0
                        j2 = 1
                    }
                } else {
                    if (i and 2 != 0) {
                        i1 = 2
                        j1 = 1
                    } else {
                        i1 = 0
                        j1 = 1
                    }
                    if (i and 4 != 0) {
                        i2 = 1
                        j2 = 2
                    } else {
                        i2 = 1
                        j2 = 0
                    }
                }
                LOOKUP_2D[i * 4 + 0] = LatticePoint2D(0, 0)
                LOOKUP_2D[i * 4 + 1] = LatticePoint2D(1, 1)
                LOOKUP_2D[i * 4 + 2] = LatticePoint2D(i1, j1)
                LOOKUP_2D[i * 4 + 3] = LatticePoint2D(i2, j2)
            }
            for (i in 0..7) {
                var i1: Int
                var j1: Int
                var k1: Int
                var i2: Int
                var j2: Int
                var k2: Int
                i1 = i shr 0 and 1
                j1 = i shr 1 and 1
                k1 = i shr 2 and 1
                i2 = i1 xor 1
                j2 = j1 xor 1
                k2 = k1 xor 1

                // The two points within this octant, one from each of the two cubic half-lattices.
                val c0 = LatticePoint3D(i1, j1, k1, 0)
                val c1 = LatticePoint3D(i1 + i2, j1 + j2, k1 + k2, 1)

                // (1, 0, 0) vs (0, 1, 1) away from octant.
                val c2 = LatticePoint3D(i1 xor 1, j1, k1, 0)
                val c3 = LatticePoint3D(i1, j1 xor 1, k1 xor 1, 0)

                // (1, 0, 0) vs (0, 1, 1) away from octant, on second half-lattice.
                val c4 = LatticePoint3D(i1 + (i2 xor 1), j1 + j2, k1 + k2, 1)
                val c5 = LatticePoint3D(i1 + i2, j1 + (j2 xor 1), k1 + (k2 xor 1), 1)

                // (0, 1, 0) vs (1, 0, 1) away from octant.
                val c6 = LatticePoint3D(i1, j1 xor 1, k1, 0)
                val c7 = LatticePoint3D(i1 xor 1, j1, k1 xor 1, 0)

                // (0, 1, 0) vs (1, 0, 1) away from octant, on second half-lattice.
                val c8 = LatticePoint3D(i1 + i2, j1 + (j2 xor 1), k1 + k2, 1)
                val c9 = LatticePoint3D(i1 + (i2 xor 1), j1 + j2, k1 + (k2 xor 1), 1)

                // (0, 0, 1) vs (1, 1, 0) away from octant.
                val cA = LatticePoint3D(i1, j1, k1 xor 1, 0)
                val cB = LatticePoint3D(i1 xor 1, j1 xor 1, k1, 0)

                // (0, 0, 1) vs (1, 1, 0) away from octant, on second half-lattice.
                val cC = LatticePoint3D(i1 + i2, j1 + j2, k1 + (k2 xor 1), 1)
                val cD = LatticePoint3D(i1 + (i2 xor 1), j1 + (j2 xor 1), k1 + k2, 1)

                // First two points are guaranteed.
                c0.nextOnSuccess = c1
                c0.nextOnFailure = c0.nextOnSuccess
                c1.nextOnSuccess = c2
                c1.nextOnFailure = c1.nextOnSuccess

                // If c2 is in range, then we know c3 and c4 are not.
                c2.nextOnFailure = c3
                c2.nextOnSuccess = c5
                c3.nextOnFailure = c4
                c3.nextOnSuccess = c4

                // If c4 is in range, then we know c5 is not.
                c4.nextOnFailure = c5
                c4.nextOnSuccess = c6
                c5.nextOnSuccess = c6
                c5.nextOnFailure = c5.nextOnSuccess

                // If c6 is in range, then we know c7 and c8 are not.
                c6.nextOnFailure = c7
                c6.nextOnSuccess = c9
                c7.nextOnFailure = c8
                c7.nextOnSuccess = c8

                // If c8 is in range, then we know c9 is not.
                c8.nextOnFailure = c9
                c8.nextOnSuccess = cA
                c9.nextOnSuccess = cA
                c9.nextOnFailure = c9.nextOnSuccess

                // If cA is in range, then we know cB and cC are not.
                cA.nextOnFailure = cB
                cA.nextOnSuccess = cD
                cB.nextOnFailure = cC
                cB.nextOnSuccess = cC

                // If cC is in range, then we know cD is not.
                cC.nextOnFailure = cD
                cC.nextOnSuccess = null
                cD.nextOnSuccess = null
                cD.nextOnFailure = cD.nextOnSuccess
                LOOKUP_3D[i] = c0
            }
            val lookup4DPregen = arrayOf(
                intArrayOf(
                    0x15,
                    0x45,
                    0x51,
                    0x54,
                    0x55,
                    0x56,
                    0x59,
                    0x5A,
                    0x65,
                    0x66,
                    0x69,
                    0x6A,
                    0x95,
                    0x96,
                    0x99,
                    0x9A,
                    0xA5,
                    0xA6,
                    0xA9,
                    0xAA
                ),
                intArrayOf(0x15, 0x45, 0x51, 0x55, 0x56, 0x59, 0x5A, 0x65, 0x66, 0x6A, 0x95, 0x96, 0x9A, 0xA6, 0xAA),
                intArrayOf(
                    0x01,
                    0x05,
                    0x11,
                    0x15,
                    0x41,
                    0x45,
                    0x51,
                    0x55,
                    0x56,
                    0x5A,
                    0x66,
                    0x6A,
                    0x96,
                    0x9A,
                    0xA6,
                    0xAA
                ),
                intArrayOf(
                    0x01,
                    0x15,
                    0x16,
                    0x45,
                    0x46,
                    0x51,
                    0x52,
                    0x55,
                    0x56,
                    0x5A,
                    0x66,
                    0x6A,
                    0x96,
                    0x9A,
                    0xA6,
                    0xAA,
                    0xAB
                ),
                intArrayOf(0x15, 0x45, 0x54, 0x55, 0x56, 0x59, 0x5A, 0x65, 0x69, 0x6A, 0x95, 0x99, 0x9A, 0xA9, 0xAA),
                intArrayOf(
                    0x05,
                    0x15,
                    0x45,
                    0x55,
                    0x56,
                    0x59,
                    0x5A,
                    0x65,
                    0x66,
                    0x69,
                    0x6A,
                    0x95,
                    0x96,
                    0x99,
                    0x9A,
                    0xAA
                ),
                intArrayOf(0x05, 0x15, 0x45, 0x55, 0x56, 0x59, 0x5A, 0x66, 0x6A, 0x96, 0x9A, 0xAA),
                intArrayOf(0x05, 0x15, 0x16, 0x45, 0x46, 0x55, 0x56, 0x59, 0x5A, 0x66, 0x6A, 0x96, 0x9A, 0xAA, 0xAB),
                intArrayOf(
                    0x04,
                    0x05,
                    0x14,
                    0x15,
                    0x44,
                    0x45,
                    0x54,
                    0x55,
                    0x59,
                    0x5A,
                    0x69,
                    0x6A,
                    0x99,
                    0x9A,
                    0xA9,
                    0xAA
                ),
                intArrayOf(0x05, 0x15, 0x45, 0x55, 0x56, 0x59, 0x5A, 0x69, 0x6A, 0x99, 0x9A, 0xAA),
                intArrayOf(0x05, 0x15, 0x45, 0x55, 0x56, 0x59, 0x5A, 0x6A, 0x9A, 0xAA),
                intArrayOf(0x05, 0x15, 0x16, 0x45, 0x46, 0x55, 0x56, 0x59, 0x5A, 0x5B, 0x6A, 0x9A, 0xAA, 0xAB),
                intArrayOf(
                    0x04,
                    0x15,
                    0x19,
                    0x45,
                    0x49,
                    0x54,
                    0x55,
                    0x58,
                    0x59,
                    0x5A,
                    0x69,
                    0x6A,
                    0x99,
                    0x9A,
                    0xA9,
                    0xAA,
                    0xAE
                ),
                intArrayOf(0x05, 0x15, 0x19, 0x45, 0x49, 0x55, 0x56, 0x59, 0x5A, 0x69, 0x6A, 0x99, 0x9A, 0xAA, 0xAE),
                intArrayOf(0x05, 0x15, 0x19, 0x45, 0x49, 0x55, 0x56, 0x59, 0x5A, 0x5E, 0x6A, 0x9A, 0xAA, 0xAE),
                intArrayOf(
                    0x05,
                    0x15,
                    0x1A,
                    0x45,
                    0x4A,
                    0x55,
                    0x56,
                    0x59,
                    0x5A,
                    0x5B,
                    0x5E,
                    0x6A,
                    0x9A,
                    0xAA,
                    0xAB,
                    0xAE,
                    0xAF
                ),
                intArrayOf(0x15, 0x51, 0x54, 0x55, 0x56, 0x59, 0x65, 0x66, 0x69, 0x6A, 0x95, 0xA5, 0xA6, 0xA9, 0xAA),
                intArrayOf(
                    0x11,
                    0x15,
                    0x51,
                    0x55,
                    0x56,
                    0x59,
                    0x5A,
                    0x65,
                    0x66,
                    0x69,
                    0x6A,
                    0x95,
                    0x96,
                    0xA5,
                    0xA6,
                    0xAA
                ),
                intArrayOf(0x11, 0x15, 0x51, 0x55, 0x56, 0x5A, 0x65, 0x66, 0x6A, 0x96, 0xA6, 0xAA),
                intArrayOf(0x11, 0x15, 0x16, 0x51, 0x52, 0x55, 0x56, 0x5A, 0x65, 0x66, 0x6A, 0x96, 0xA6, 0xAA, 0xAB),
                intArrayOf(
                    0x14,
                    0x15,
                    0x54,
                    0x55,
                    0x56,
                    0x59,
                    0x5A,
                    0x65,
                    0x66,
                    0x69,
                    0x6A,
                    0x95,
                    0x99,
                    0xA5,
                    0xA9,
                    0xAA
                ),
                intArrayOf(0x15, 0x55, 0x56, 0x59, 0x5A, 0x65, 0x66, 0x69, 0x6A, 0x95, 0x9A, 0xA6, 0xA9, 0xAA),
                intArrayOf(0x15, 0x55, 0x56, 0x59, 0x5A, 0x65, 0x66, 0x69, 0x6A, 0x96, 0x9A, 0xA6, 0xAA, 0xAB),
                intArrayOf(0x15, 0x16, 0x55, 0x56, 0x5A, 0x66, 0x6A, 0x6B, 0x96, 0x9A, 0xA6, 0xAA, 0xAB),
                intArrayOf(0x14, 0x15, 0x54, 0x55, 0x59, 0x5A, 0x65, 0x69, 0x6A, 0x99, 0xA9, 0xAA),
                intArrayOf(0x15, 0x55, 0x56, 0x59, 0x5A, 0x65, 0x66, 0x69, 0x6A, 0x99, 0x9A, 0xA9, 0xAA, 0xAE),
                intArrayOf(0x15, 0x55, 0x56, 0x59, 0x5A, 0x65, 0x66, 0x69, 0x6A, 0x9A, 0xAA),
                intArrayOf(0x15, 0x16, 0x55, 0x56, 0x59, 0x5A, 0x66, 0x6A, 0x6B, 0x9A, 0xAA, 0xAB),
                intArrayOf(0x14, 0x15, 0x19, 0x54, 0x55, 0x58, 0x59, 0x5A, 0x65, 0x69, 0x6A, 0x99, 0xA9, 0xAA, 0xAE),
                intArrayOf(0x15, 0x19, 0x55, 0x59, 0x5A, 0x69, 0x6A, 0x6E, 0x99, 0x9A, 0xA9, 0xAA, 0xAE),
                intArrayOf(0x15, 0x19, 0x55, 0x56, 0x59, 0x5A, 0x69, 0x6A, 0x6E, 0x9A, 0xAA, 0xAE),
                intArrayOf(0x15, 0x1A, 0x55, 0x56, 0x59, 0x5A, 0x6A, 0x6B, 0x6E, 0x9A, 0xAA, 0xAB, 0xAE, 0xAF),
                intArrayOf(
                    0x10,
                    0x11,
                    0x14,
                    0x15,
                    0x50,
                    0x51,
                    0x54,
                    0x55,
                    0x65,
                    0x66,
                    0x69,
                    0x6A,
                    0xA5,
                    0xA6,
                    0xA9,
                    0xAA
                ),
                intArrayOf(0x11, 0x15, 0x51, 0x55, 0x56, 0x65, 0x66, 0x69, 0x6A, 0xA5, 0xA6, 0xAA),
                intArrayOf(0x11, 0x15, 0x51, 0x55, 0x56, 0x65, 0x66, 0x6A, 0xA6, 0xAA),
                intArrayOf(0x11, 0x15, 0x16, 0x51, 0x52, 0x55, 0x56, 0x65, 0x66, 0x67, 0x6A, 0xA6, 0xAA, 0xAB),
                intArrayOf(0x14, 0x15, 0x54, 0x55, 0x59, 0x65, 0x66, 0x69, 0x6A, 0xA5, 0xA9, 0xAA),
                intArrayOf(0x15, 0x55, 0x56, 0x59, 0x5A, 0x65, 0x66, 0x69, 0x6A, 0xA5, 0xA6, 0xA9, 0xAA, 0xBA),
                intArrayOf(0x15, 0x55, 0x56, 0x59, 0x5A, 0x65, 0x66, 0x69, 0x6A, 0xA6, 0xAA),
                intArrayOf(0x15, 0x16, 0x55, 0x56, 0x5A, 0x65, 0x66, 0x6A, 0x6B, 0xA6, 0xAA, 0xAB),
                intArrayOf(0x14, 0x15, 0x54, 0x55, 0x59, 0x65, 0x69, 0x6A, 0xA9, 0xAA),
                intArrayOf(0x15, 0x55, 0x56, 0x59, 0x5A, 0x65, 0x66, 0x69, 0x6A, 0xA9, 0xAA),
                intArrayOf(0x15, 0x55, 0x56, 0x59, 0x5A, 0x65, 0x66, 0x69, 0x6A, 0xAA),
                intArrayOf(0x15, 0x16, 0x55, 0x56, 0x59, 0x5A, 0x65, 0x66, 0x69, 0x6A, 0x6B, 0xAA, 0xAB),
                intArrayOf(0x14, 0x15, 0x19, 0x54, 0x55, 0x58, 0x59, 0x65, 0x69, 0x6A, 0x6D, 0xA9, 0xAA, 0xAE),
                intArrayOf(0x15, 0x19, 0x55, 0x59, 0x5A, 0x65, 0x69, 0x6A, 0x6E, 0xA9, 0xAA, 0xAE),
                intArrayOf(0x15, 0x19, 0x55, 0x56, 0x59, 0x5A, 0x65, 0x66, 0x69, 0x6A, 0x6E, 0xAA, 0xAE),
                intArrayOf(0x15, 0x55, 0x56, 0x59, 0x5A, 0x66, 0x69, 0x6A, 0x6B, 0x6E, 0x9A, 0xAA, 0xAB, 0xAE, 0xAF),
                intArrayOf(
                    0x10,
                    0x15,
                    0x25,
                    0x51,
                    0x54,
                    0x55,
                    0x61,
                    0x64,
                    0x65,
                    0x66,
                    0x69,
                    0x6A,
                    0xA5,
                    0xA6,
                    0xA9,
                    0xAA,
                    0xBA
                ),
                intArrayOf(0x11, 0x15, 0x25, 0x51, 0x55, 0x56, 0x61, 0x65, 0x66, 0x69, 0x6A, 0xA5, 0xA6, 0xAA, 0xBA),
                intArrayOf(0x11, 0x15, 0x25, 0x51, 0x55, 0x56, 0x61, 0x65, 0x66, 0x6A, 0x76, 0xA6, 0xAA, 0xBA),
                intArrayOf(
                    0x11,
                    0x15,
                    0x26,
                    0x51,
                    0x55,
                    0x56,
                    0x62,
                    0x65,
                    0x66,
                    0x67,
                    0x6A,
                    0x76,
                    0xA6,
                    0xAA,
                    0xAB,
                    0xBA,
                    0xBB
                ),
                intArrayOf(0x14, 0x15, 0x25, 0x54, 0x55, 0x59, 0x64, 0x65, 0x66, 0x69, 0x6A, 0xA5, 0xA9, 0xAA, 0xBA),
                intArrayOf(0x15, 0x25, 0x55, 0x65, 0x66, 0x69, 0x6A, 0x7A, 0xA5, 0xA6, 0xA9, 0xAA, 0xBA),
                intArrayOf(0x15, 0x25, 0x55, 0x56, 0x65, 0x66, 0x69, 0x6A, 0x7A, 0xA6, 0xAA, 0xBA),
                intArrayOf(0x15, 0x26, 0x55, 0x56, 0x65, 0x66, 0x6A, 0x6B, 0x7A, 0xA6, 0xAA, 0xAB, 0xBA, 0xBB),
                intArrayOf(0x14, 0x15, 0x25, 0x54, 0x55, 0x59, 0x64, 0x65, 0x69, 0x6A, 0x79, 0xA9, 0xAA, 0xBA),
                intArrayOf(0x15, 0x25, 0x55, 0x59, 0x65, 0x66, 0x69, 0x6A, 0x7A, 0xA9, 0xAA, 0xBA),
                intArrayOf(0x15, 0x25, 0x55, 0x56, 0x59, 0x5A, 0x65, 0x66, 0x69, 0x6A, 0x7A, 0xAA, 0xBA),
                intArrayOf(0x15, 0x55, 0x56, 0x5A, 0x65, 0x66, 0x69, 0x6A, 0x6B, 0x7A, 0xA6, 0xAA, 0xAB, 0xBA, 0xBB),
                intArrayOf(
                    0x14,
                    0x15,
                    0x29,
                    0x54,
                    0x55,
                    0x59,
                    0x65,
                    0x68,
                    0x69,
                    0x6A,
                    0x6D,
                    0x79,
                    0xA9,
                    0xAA,
                    0xAE,
                    0xBA,
                    0xBE
                ),
                intArrayOf(0x15, 0x29, 0x55, 0x59, 0x65, 0x69, 0x6A, 0x6E, 0x7A, 0xA9, 0xAA, 0xAE, 0xBA, 0xBE),
                intArrayOf(0x15, 0x55, 0x59, 0x5A, 0x65, 0x66, 0x69, 0x6A, 0x6E, 0x7A, 0xA9, 0xAA, 0xAE, 0xBA, 0xBE),
                intArrayOf(
                    0x15,
                    0x55,
                    0x56,
                    0x59,
                    0x5A,
                    0x65,
                    0x66,
                    0x69,
                    0x6A,
                    0x6B,
                    0x6E,
                    0x7A,
                    0xAA,
                    0xAB,
                    0xAE,
                    0xBA,
                    0xBF
                ),
                intArrayOf(0x45, 0x51, 0x54, 0x55, 0x56, 0x59, 0x65, 0x95, 0x96, 0x99, 0x9A, 0xA5, 0xA6, 0xA9, 0xAA),
                intArrayOf(
                    0x41,
                    0x45,
                    0x51,
                    0x55,
                    0x56,
                    0x59,
                    0x5A,
                    0x65,
                    0x66,
                    0x95,
                    0x96,
                    0x99,
                    0x9A,
                    0xA5,
                    0xA6,
                    0xAA
                ),
                intArrayOf(0x41, 0x45, 0x51, 0x55, 0x56, 0x5A, 0x66, 0x95, 0x96, 0x9A, 0xA6, 0xAA),
                intArrayOf(0x41, 0x45, 0x46, 0x51, 0x52, 0x55, 0x56, 0x5A, 0x66, 0x95, 0x96, 0x9A, 0xA6, 0xAA, 0xAB),
                intArrayOf(
                    0x44,
                    0x45,
                    0x54,
                    0x55,
                    0x56,
                    0x59,
                    0x5A,
                    0x65,
                    0x69,
                    0x95,
                    0x96,
                    0x99,
                    0x9A,
                    0xA5,
                    0xA9,
                    0xAA
                ),
                intArrayOf(0x45, 0x55, 0x56, 0x59, 0x5A, 0x65, 0x6A, 0x95, 0x96, 0x99, 0x9A, 0xA6, 0xA9, 0xAA),
                intArrayOf(0x45, 0x55, 0x56, 0x59, 0x5A, 0x66, 0x6A, 0x95, 0x96, 0x99, 0x9A, 0xA6, 0xAA, 0xAB),
                intArrayOf(0x45, 0x46, 0x55, 0x56, 0x5A, 0x66, 0x6A, 0x96, 0x9A, 0x9B, 0xA6, 0xAA, 0xAB),
                intArrayOf(0x44, 0x45, 0x54, 0x55, 0x59, 0x5A, 0x69, 0x95, 0x99, 0x9A, 0xA9, 0xAA),
                intArrayOf(0x45, 0x55, 0x56, 0x59, 0x5A, 0x69, 0x6A, 0x95, 0x96, 0x99, 0x9A, 0xA9, 0xAA, 0xAE),
                intArrayOf(0x45, 0x55, 0x56, 0x59, 0x5A, 0x6A, 0x95, 0x96, 0x99, 0x9A, 0xAA),
                intArrayOf(0x45, 0x46, 0x55, 0x56, 0x59, 0x5A, 0x6A, 0x96, 0x9A, 0x9B, 0xAA, 0xAB),
                intArrayOf(0x44, 0x45, 0x49, 0x54, 0x55, 0x58, 0x59, 0x5A, 0x69, 0x95, 0x99, 0x9A, 0xA9, 0xAA, 0xAE),
                intArrayOf(0x45, 0x49, 0x55, 0x59, 0x5A, 0x69, 0x6A, 0x99, 0x9A, 0x9E, 0xA9, 0xAA, 0xAE),
                intArrayOf(0x45, 0x49, 0x55, 0x56, 0x59, 0x5A, 0x6A, 0x99, 0x9A, 0x9E, 0xAA, 0xAE),
                intArrayOf(0x45, 0x4A, 0x55, 0x56, 0x59, 0x5A, 0x6A, 0x9A, 0x9B, 0x9E, 0xAA, 0xAB, 0xAE, 0xAF),
                intArrayOf(
                    0x50,
                    0x51,
                    0x54,
                    0x55,
                    0x56,
                    0x59,
                    0x65,
                    0x66,
                    0x69,
                    0x95,
                    0x96,
                    0x99,
                    0xA5,
                    0xA6,
                    0xA9,
                    0xAA
                ),
                intArrayOf(0x51, 0x55, 0x56, 0x59, 0x65, 0x66, 0x6A, 0x95, 0x96, 0x9A, 0xA5, 0xA6, 0xA9, 0xAA),
                intArrayOf(0x51, 0x55, 0x56, 0x5A, 0x65, 0x66, 0x6A, 0x95, 0x96, 0x9A, 0xA5, 0xA6, 0xAA, 0xAB),
                intArrayOf(0x51, 0x52, 0x55, 0x56, 0x5A, 0x66, 0x6A, 0x96, 0x9A, 0xA6, 0xA7, 0xAA, 0xAB),
                intArrayOf(0x54, 0x55, 0x56, 0x59, 0x65, 0x69, 0x6A, 0x95, 0x99, 0x9A, 0xA5, 0xA6, 0xA9, 0xAA),
                intArrayOf(
                    0x55,
                    0x56,
                    0x59,
                    0x5A,
                    0x65,
                    0x66,
                    0x69,
                    0x6A,
                    0x95,
                    0x96,
                    0x99,
                    0x9A,
                    0xA5,
                    0xA6,
                    0xA9,
                    0xAA
                ),
                intArrayOf(
                    0x15,
                    0x45,
                    0x51,
                    0x55,
                    0x56,
                    0x59,
                    0x5A,
                    0x65,
                    0x66,
                    0x6A,
                    0x95,
                    0x96,
                    0x9A,
                    0xA6,
                    0xAA,
                    0xAB
                ),
                intArrayOf(0x55, 0x56, 0x5A, 0x66, 0x6A, 0x96, 0x9A, 0xA6, 0xAA, 0xAB),
                intArrayOf(0x54, 0x55, 0x59, 0x5A, 0x65, 0x69, 0x6A, 0x95, 0x99, 0x9A, 0xA5, 0xA9, 0xAA, 0xAE),
                intArrayOf(
                    0x15,
                    0x45,
                    0x54,
                    0x55,
                    0x56,
                    0x59,
                    0x5A,
                    0x65,
                    0x69,
                    0x6A,
                    0x95,
                    0x99,
                    0x9A,
                    0xA9,
                    0xAA,
                    0xAE
                ),
                intArrayOf(
                    0x15,
                    0x45,
                    0x55,
                    0x56,
                    0x59,
                    0x5A,
                    0x65,
                    0x66,
                    0x69,
                    0x6A,
                    0x95,
                    0x96,
                    0x99,
                    0x9A,
                    0xA6,
                    0xA9,
                    0xAA,
                    0xAB,
                    0xAE
                ),
                intArrayOf(0x55, 0x56, 0x59, 0x5A, 0x66, 0x6A, 0x96, 0x9A, 0xA6, 0xAA, 0xAB),
                intArrayOf(0x54, 0x55, 0x58, 0x59, 0x5A, 0x69, 0x6A, 0x99, 0x9A, 0xA9, 0xAA, 0xAD, 0xAE),
                intArrayOf(0x55, 0x59, 0x5A, 0x69, 0x6A, 0x99, 0x9A, 0xA9, 0xAA, 0xAE),
                intArrayOf(0x55, 0x56, 0x59, 0x5A, 0x69, 0x6A, 0x99, 0x9A, 0xA9, 0xAA, 0xAE),
                intArrayOf(0x55, 0x56, 0x59, 0x5A, 0x6A, 0x9A, 0xAA, 0xAB, 0xAE, 0xAF),
                intArrayOf(0x50, 0x51, 0x54, 0x55, 0x65, 0x66, 0x69, 0x95, 0xA5, 0xA6, 0xA9, 0xAA),
                intArrayOf(0x51, 0x55, 0x56, 0x65, 0x66, 0x69, 0x6A, 0x95, 0x96, 0xA5, 0xA6, 0xA9, 0xAA, 0xBA),
                intArrayOf(0x51, 0x55, 0x56, 0x65, 0x66, 0x6A, 0x95, 0x96, 0xA5, 0xA6, 0xAA),
                intArrayOf(0x51, 0x52, 0x55, 0x56, 0x65, 0x66, 0x6A, 0x96, 0xA6, 0xA7, 0xAA, 0xAB),
                intArrayOf(0x54, 0x55, 0x59, 0x65, 0x66, 0x69, 0x6A, 0x95, 0x99, 0xA5, 0xA6, 0xA9, 0xAA, 0xBA),
                intArrayOf(
                    0x15,
                    0x51,
                    0x54,
                    0x55,
                    0x56,
                    0x59,
                    0x65,
                    0x66,
                    0x69,
                    0x6A,
                    0x95,
                    0xA5,
                    0xA6,
                    0xA9,
                    0xAA,
                    0xBA
                ),
                intArrayOf(
                    0x15,
                    0x51,
                    0x55,
                    0x56,
                    0x59,
                    0x5A,
                    0x65,
                    0x66,
                    0x69,
                    0x6A,
                    0x95,
                    0x96,
                    0x9A,
                    0xA5,
                    0xA6,
                    0xA9,
                    0xAA,
                    0xAB,
                    0xBA
                ),
                intArrayOf(0x55, 0x56, 0x5A, 0x65, 0x66, 0x6A, 0x96, 0x9A, 0xA6, 0xAA, 0xAB),
                intArrayOf(0x54, 0x55, 0x59, 0x65, 0x69, 0x6A, 0x95, 0x99, 0xA5, 0xA9, 0xAA),
                intArrayOf(
                    0x15,
                    0x54,
                    0x55,
                    0x56,
                    0x59,
                    0x5A,
                    0x65,
                    0x66,
                    0x69,
                    0x6A,
                    0x95,
                    0x99,
                    0x9A,
                    0xA5,
                    0xA6,
                    0xA9,
                    0xAA,
                    0xAE,
                    0xBA
                ),
                intArrayOf(0x15, 0x55, 0x56, 0x59, 0x5A, 0x65, 0x66, 0x69, 0x6A, 0x9A, 0xA6, 0xA9, 0xAA),
                intArrayOf(0x15, 0x55, 0x56, 0x59, 0x5A, 0x65, 0x66, 0x69, 0x6A, 0x96, 0x9A, 0xA6, 0xAA, 0xAB),
                intArrayOf(0x54, 0x55, 0x58, 0x59, 0x65, 0x69, 0x6A, 0x99, 0xA9, 0xAA, 0xAD, 0xAE),
                intArrayOf(0x55, 0x59, 0x5A, 0x65, 0x69, 0x6A, 0x99, 0x9A, 0xA9, 0xAA, 0xAE),
                intArrayOf(0x15, 0x55, 0x56, 0x59, 0x5A, 0x65, 0x66, 0x69, 0x6A, 0x99, 0x9A, 0xA9, 0xAA, 0xAE),
                intArrayOf(0x15, 0x55, 0x56, 0x59, 0x5A, 0x66, 0x69, 0x6A, 0x9A, 0xAA, 0xAB, 0xAE, 0xAF),
                intArrayOf(0x50, 0x51, 0x54, 0x55, 0x61, 0x64, 0x65, 0x66, 0x69, 0x95, 0xA5, 0xA6, 0xA9, 0xAA, 0xBA),
                intArrayOf(0x51, 0x55, 0x61, 0x65, 0x66, 0x69, 0x6A, 0xA5, 0xA6, 0xA9, 0xAA, 0xB6, 0xBA),
                intArrayOf(0x51, 0x55, 0x56, 0x61, 0x65, 0x66, 0x6A, 0xA5, 0xA6, 0xAA, 0xB6, 0xBA),
                intArrayOf(0x51, 0x55, 0x56, 0x62, 0x65, 0x66, 0x6A, 0xA6, 0xA7, 0xAA, 0xAB, 0xB6, 0xBA, 0xBB),
                intArrayOf(0x54, 0x55, 0x64, 0x65, 0x66, 0x69, 0x6A, 0xA5, 0xA6, 0xA9, 0xAA, 0xB9, 0xBA),
                intArrayOf(0x55, 0x65, 0x66, 0x69, 0x6A, 0xA5, 0xA6, 0xA9, 0xAA, 0xBA),
                intArrayOf(0x55, 0x56, 0x65, 0x66, 0x69, 0x6A, 0xA5, 0xA6, 0xA9, 0xAA, 0xBA),
                intArrayOf(0x55, 0x56, 0x65, 0x66, 0x6A, 0xA6, 0xAA, 0xAB, 0xBA, 0xBB),
                intArrayOf(0x54, 0x55, 0x59, 0x64, 0x65, 0x69, 0x6A, 0xA5, 0xA9, 0xAA, 0xB9, 0xBA),
                intArrayOf(0x55, 0x59, 0x65, 0x66, 0x69, 0x6A, 0xA5, 0xA6, 0xA9, 0xAA, 0xBA),
                intArrayOf(0x15, 0x55, 0x56, 0x59, 0x5A, 0x65, 0x66, 0x69, 0x6A, 0xA5, 0xA6, 0xA9, 0xAA, 0xBA),
                intArrayOf(0x15, 0x55, 0x56, 0x5A, 0x65, 0x66, 0x69, 0x6A, 0xA6, 0xAA, 0xAB, 0xBA, 0xBB),
                intArrayOf(0x54, 0x55, 0x59, 0x65, 0x68, 0x69, 0x6A, 0xA9, 0xAA, 0xAD, 0xAE, 0xB9, 0xBA, 0xBE),
                intArrayOf(0x55, 0x59, 0x65, 0x69, 0x6A, 0xA9, 0xAA, 0xAE, 0xBA, 0xBE),
                intArrayOf(0x15, 0x55, 0x59, 0x5A, 0x65, 0x66, 0x69, 0x6A, 0xA9, 0xAA, 0xAE, 0xBA, 0xBE),
                intArrayOf(0x55, 0x56, 0x59, 0x5A, 0x65, 0x66, 0x69, 0x6A, 0xAA, 0xAB, 0xAE, 0xBA, 0xBF),
                intArrayOf(
                    0x40,
                    0x41,
                    0x44,
                    0x45,
                    0x50,
                    0x51,
                    0x54,
                    0x55,
                    0x95,
                    0x96,
                    0x99,
                    0x9A,
                    0xA5,
                    0xA6,
                    0xA9,
                    0xAA
                ),
                intArrayOf(0x41, 0x45, 0x51, 0x55, 0x56, 0x95, 0x96, 0x99, 0x9A, 0xA5, 0xA6, 0xAA),
                intArrayOf(0x41, 0x45, 0x51, 0x55, 0x56, 0x95, 0x96, 0x9A, 0xA6, 0xAA),
                intArrayOf(0x41, 0x45, 0x46, 0x51, 0x52, 0x55, 0x56, 0x95, 0x96, 0x97, 0x9A, 0xA6, 0xAA, 0xAB),
                intArrayOf(0x44, 0x45, 0x54, 0x55, 0x59, 0x95, 0x96, 0x99, 0x9A, 0xA5, 0xA9, 0xAA),
                intArrayOf(0x45, 0x55, 0x56, 0x59, 0x5A, 0x95, 0x96, 0x99, 0x9A, 0xA5, 0xA6, 0xA9, 0xAA, 0xEA),
                intArrayOf(0x45, 0x55, 0x56, 0x59, 0x5A, 0x95, 0x96, 0x99, 0x9A, 0xA6, 0xAA),
                intArrayOf(0x45, 0x46, 0x55, 0x56, 0x5A, 0x95, 0x96, 0x9A, 0x9B, 0xA6, 0xAA, 0xAB),
                intArrayOf(0x44, 0x45, 0x54, 0x55, 0x59, 0x95, 0x99, 0x9A, 0xA9, 0xAA),
                intArrayOf(0x45, 0x55, 0x56, 0x59, 0x5A, 0x95, 0x96, 0x99, 0x9A, 0xA9, 0xAA),
                intArrayOf(0x45, 0x55, 0x56, 0x59, 0x5A, 0x95, 0x96, 0x99, 0x9A, 0xAA),
                intArrayOf(0x45, 0x46, 0x55, 0x56, 0x59, 0x5A, 0x95, 0x96, 0x99, 0x9A, 0x9B, 0xAA, 0xAB),
                intArrayOf(0x44, 0x45, 0x49, 0x54, 0x55, 0x58, 0x59, 0x95, 0x99, 0x9A, 0x9D, 0xA9, 0xAA, 0xAE),
                intArrayOf(0x45, 0x49, 0x55, 0x59, 0x5A, 0x95, 0x99, 0x9A, 0x9E, 0xA9, 0xAA, 0xAE),
                intArrayOf(0x45, 0x49, 0x55, 0x56, 0x59, 0x5A, 0x95, 0x96, 0x99, 0x9A, 0x9E, 0xAA, 0xAE),
                intArrayOf(0x45, 0x55, 0x56, 0x59, 0x5A, 0x6A, 0x96, 0x99, 0x9A, 0x9B, 0x9E, 0xAA, 0xAB, 0xAE, 0xAF),
                intArrayOf(0x50, 0x51, 0x54, 0x55, 0x65, 0x95, 0x96, 0x99, 0xA5, 0xA6, 0xA9, 0xAA),
                intArrayOf(0x51, 0x55, 0x56, 0x65, 0x66, 0x95, 0x96, 0x99, 0x9A, 0xA5, 0xA6, 0xA9, 0xAA, 0xEA),
                intArrayOf(0x51, 0x55, 0x56, 0x65, 0x66, 0x95, 0x96, 0x9A, 0xA5, 0xA6, 0xAA),
                intArrayOf(0x51, 0x52, 0x55, 0x56, 0x66, 0x95, 0x96, 0x9A, 0xA6, 0xA7, 0xAA, 0xAB),
                intArrayOf(0x54, 0x55, 0x59, 0x65, 0x69, 0x95, 0x96, 0x99, 0x9A, 0xA5, 0xA6, 0xA9, 0xAA, 0xEA),
                intArrayOf(
                    0x45,
                    0x51,
                    0x54,
                    0x55,
                    0x56,
                    0x59,
                    0x65,
                    0x95,
                    0x96,
                    0x99,
                    0x9A,
                    0xA5,
                    0xA6,
                    0xA9,
                    0xAA,
                    0xEA
                ),
                intArrayOf(
                    0x45,
                    0x51,
                    0x55,
                    0x56,
                    0x59,
                    0x5A,
                    0x65,
                    0x66,
                    0x6A,
                    0x95,
                    0x96,
                    0x99,
                    0x9A,
                    0xA5,
                    0xA6,
                    0xA9,
                    0xAA,
                    0xAB,
                    0xEA
                ),
                intArrayOf(0x55, 0x56, 0x5A, 0x66, 0x6A, 0x95, 0x96, 0x9A, 0xA6, 0xAA, 0xAB),
                intArrayOf(0x54, 0x55, 0x59, 0x65, 0x69, 0x95, 0x99, 0x9A, 0xA5, 0xA9, 0xAA),
                intArrayOf(
                    0x45,
                    0x54,
                    0x55,
                    0x56,
                    0x59,
                    0x5A,
                    0x65,
                    0x69,
                    0x6A,
                    0x95,
                    0x96,
                    0x99,
                    0x9A,
                    0xA5,
                    0xA6,
                    0xA9,
                    0xAA,
                    0xAE,
                    0xEA
                ),
                intArrayOf(0x45, 0x55, 0x56, 0x59, 0x5A, 0x6A, 0x95, 0x96, 0x99, 0x9A, 0xA6, 0xA9, 0xAA),
                intArrayOf(0x45, 0x55, 0x56, 0x59, 0x5A, 0x66, 0x6A, 0x95, 0x96, 0x99, 0x9A, 0xA6, 0xAA, 0xAB),
                intArrayOf(0x54, 0x55, 0x58, 0x59, 0x69, 0x95, 0x99, 0x9A, 0xA9, 0xAA, 0xAD, 0xAE),
                intArrayOf(0x55, 0x59, 0x5A, 0x69, 0x6A, 0x95, 0x99, 0x9A, 0xA9, 0xAA, 0xAE),
                intArrayOf(0x45, 0x55, 0x56, 0x59, 0x5A, 0x69, 0x6A, 0x95, 0x96, 0x99, 0x9A, 0xA9, 0xAA, 0xAE),
                intArrayOf(0x45, 0x55, 0x56, 0x59, 0x5A, 0x6A, 0x96, 0x99, 0x9A, 0xAA, 0xAB, 0xAE, 0xAF),
                intArrayOf(0x50, 0x51, 0x54, 0x55, 0x65, 0x95, 0xA5, 0xA6, 0xA9, 0xAA),
                intArrayOf(0x51, 0x55, 0x56, 0x65, 0x66, 0x95, 0x96, 0xA5, 0xA6, 0xA9, 0xAA),
                intArrayOf(0x51, 0x55, 0x56, 0x65, 0x66, 0x95, 0x96, 0xA5, 0xA6, 0xAA),
                intArrayOf(0x51, 0x52, 0x55, 0x56, 0x65, 0x66, 0x95, 0x96, 0xA5, 0xA6, 0xA7, 0xAA, 0xAB),
                intArrayOf(0x54, 0x55, 0x59, 0x65, 0x69, 0x95, 0x99, 0xA5, 0xA6, 0xA9, 0xAA),
                intArrayOf(
                    0x51,
                    0x54,
                    0x55,
                    0x56,
                    0x59,
                    0x65,
                    0x66,
                    0x69,
                    0x6A,
                    0x95,
                    0x96,
                    0x99,
                    0x9A,
                    0xA5,
                    0xA6,
                    0xA9,
                    0xAA,
                    0xBA,
                    0xEA
                ),
                intArrayOf(0x51, 0x55, 0x56, 0x65, 0x66, 0x6A, 0x95, 0x96, 0x9A, 0xA5, 0xA6, 0xA9, 0xAA),
                intArrayOf(0x51, 0x55, 0x56, 0x5A, 0x65, 0x66, 0x6A, 0x95, 0x96, 0x9A, 0xA5, 0xA6, 0xAA, 0xAB),
                intArrayOf(0x54, 0x55, 0x59, 0x65, 0x69, 0x95, 0x99, 0xA5, 0xA9, 0xAA),
                intArrayOf(0x54, 0x55, 0x59, 0x65, 0x69, 0x6A, 0x95, 0x99, 0x9A, 0xA5, 0xA6, 0xA9, 0xAA),
                intArrayOf(
                    0x55,
                    0x56,
                    0x59,
                    0x5A,
                    0x65,
                    0x66,
                    0x69,
                    0x6A,
                    0x95,
                    0x96,
                    0x99,
                    0x9A,
                    0xA5,
                    0xA6,
                    0xA9,
                    0xAA
                ),
                intArrayOf(0x55, 0x56, 0x59, 0x5A, 0x65, 0x66, 0x6A, 0x95, 0x96, 0x9A, 0xA6, 0xA9, 0xAA, 0xAB),
                intArrayOf(0x54, 0x55, 0x58, 0x59, 0x65, 0x69, 0x95, 0x99, 0xA5, 0xA9, 0xAA, 0xAD, 0xAE),
                intArrayOf(0x54, 0x55, 0x59, 0x5A, 0x65, 0x69, 0x6A, 0x95, 0x99, 0x9A, 0xA5, 0xA9, 0xAA, 0xAE),
                intArrayOf(0x55, 0x56, 0x59, 0x5A, 0x65, 0x69, 0x6A, 0x95, 0x99, 0x9A, 0xA6, 0xA9, 0xAA, 0xAE),
                intArrayOf(
                    0x55,
                    0x56,
                    0x59,
                    0x5A,
                    0x66,
                    0x69,
                    0x6A,
                    0x96,
                    0x99,
                    0x9A,
                    0xA6,
                    0xA9,
                    0xAA,
                    0xAB,
                    0xAE,
                    0xAF
                ),
                intArrayOf(0x50, 0x51, 0x54, 0x55, 0x61, 0x64, 0x65, 0x95, 0xA5, 0xA6, 0xA9, 0xAA, 0xB5, 0xBA),
                intArrayOf(0x51, 0x55, 0x61, 0x65, 0x66, 0x95, 0xA5, 0xA6, 0xA9, 0xAA, 0xB6, 0xBA),
                intArrayOf(0x51, 0x55, 0x56, 0x61, 0x65, 0x66, 0x95, 0x96, 0xA5, 0xA6, 0xAA, 0xB6, 0xBA),
                intArrayOf(0x51, 0x55, 0x56, 0x65, 0x66, 0x6A, 0x96, 0xA5, 0xA6, 0xA7, 0xAA, 0xAB, 0xB6, 0xBA, 0xBB),
                intArrayOf(0x54, 0x55, 0x64, 0x65, 0x69, 0x95, 0xA5, 0xA6, 0xA9, 0xAA, 0xB9, 0xBA),
                intArrayOf(0x55, 0x65, 0x66, 0x69, 0x6A, 0x95, 0xA5, 0xA6, 0xA9, 0xAA, 0xBA),
                intArrayOf(0x51, 0x55, 0x56, 0x65, 0x66, 0x69, 0x6A, 0x95, 0x96, 0xA5, 0xA6, 0xA9, 0xAA, 0xBA),
                intArrayOf(0x51, 0x55, 0x56, 0x65, 0x66, 0x6A, 0x96, 0xA5, 0xA6, 0xAA, 0xAB, 0xBA, 0xBB),
                intArrayOf(0x54, 0x55, 0x59, 0x64, 0x65, 0x69, 0x95, 0x99, 0xA5, 0xA9, 0xAA, 0xB9, 0xBA),
                intArrayOf(0x54, 0x55, 0x59, 0x65, 0x66, 0x69, 0x6A, 0x95, 0x99, 0xA5, 0xA6, 0xA9, 0xAA, 0xBA),
                intArrayOf(0x55, 0x56, 0x59, 0x65, 0x66, 0x69, 0x6A, 0x95, 0x9A, 0xA5, 0xA6, 0xA9, 0xAA, 0xBA),
                intArrayOf(
                    0x55,
                    0x56,
                    0x5A,
                    0x65,
                    0x66,
                    0x69,
                    0x6A,
                    0x96,
                    0x9A,
                    0xA5,
                    0xA6,
                    0xA9,
                    0xAA,
                    0xAB,
                    0xBA,
                    0xBB
                ),
                intArrayOf(0x54, 0x55, 0x59, 0x65, 0x69, 0x6A, 0x99, 0xA5, 0xA9, 0xAA, 0xAD, 0xAE, 0xB9, 0xBA, 0xBE),
                intArrayOf(0x54, 0x55, 0x59, 0x65, 0x69, 0x6A, 0x99, 0xA5, 0xA9, 0xAA, 0xAE, 0xBA, 0xBE),
                intArrayOf(
                    0x55,
                    0x59,
                    0x5A,
                    0x65,
                    0x66,
                    0x69,
                    0x6A,
                    0x99,
                    0x9A,
                    0xA5,
                    0xA6,
                    0xA9,
                    0xAA,
                    0xAE,
                    0xBA,
                    0xBE
                ),
                intArrayOf(0x55, 0x56, 0x59, 0x5A, 0x65, 0x66, 0x69, 0x6A, 0x9A, 0xA6, 0xA9, 0xAA, 0xAB, 0xAE, 0xBA),
                intArrayOf(
                    0x40,
                    0x45,
                    0x51,
                    0x54,
                    0x55,
                    0x85,
                    0x91,
                    0x94,
                    0x95,
                    0x96,
                    0x99,
                    0x9A,
                    0xA5,
                    0xA6,
                    0xA9,
                    0xAA,
                    0xEA
                ),
                intArrayOf(0x41, 0x45, 0x51, 0x55, 0x56, 0x85, 0x91, 0x95, 0x96, 0x99, 0x9A, 0xA5, 0xA6, 0xAA, 0xEA),
                intArrayOf(0x41, 0x45, 0x51, 0x55, 0x56, 0x85, 0x91, 0x95, 0x96, 0x9A, 0xA6, 0xAA, 0xD6, 0xEA),
                intArrayOf(
                    0x41,
                    0x45,
                    0x51,
                    0x55,
                    0x56,
                    0x86,
                    0x92,
                    0x95,
                    0x96,
                    0x97,
                    0x9A,
                    0xA6,
                    0xAA,
                    0xAB,
                    0xD6,
                    0xEA,
                    0xEB
                ),
                intArrayOf(0x44, 0x45, 0x54, 0x55, 0x59, 0x85, 0x94, 0x95, 0x96, 0x99, 0x9A, 0xA5, 0xA9, 0xAA, 0xEA),
                intArrayOf(0x45, 0x55, 0x85, 0x95, 0x96, 0x99, 0x9A, 0xA5, 0xA6, 0xA9, 0xAA, 0xDA, 0xEA),
                intArrayOf(0x45, 0x55, 0x56, 0x85, 0x95, 0x96, 0x99, 0x9A, 0xA6, 0xAA, 0xDA, 0xEA),
                intArrayOf(0x45, 0x55, 0x56, 0x86, 0x95, 0x96, 0x9A, 0x9B, 0xA6, 0xAA, 0xAB, 0xDA, 0xEA, 0xEB),
                intArrayOf(0x44, 0x45, 0x54, 0x55, 0x59, 0x85, 0x94, 0x95, 0x99, 0x9A, 0xA9, 0xAA, 0xD9, 0xEA),
                intArrayOf(0x45, 0x55, 0x59, 0x85, 0x95, 0x96, 0x99, 0x9A, 0xA9, 0xAA, 0xDA, 0xEA),
                intArrayOf(0x45, 0x55, 0x56, 0x59, 0x5A, 0x85, 0x95, 0x96, 0x99, 0x9A, 0xAA, 0xDA, 0xEA),
                intArrayOf(0x45, 0x55, 0x56, 0x5A, 0x95, 0x96, 0x99, 0x9A, 0x9B, 0xA6, 0xAA, 0xAB, 0xDA, 0xEA, 0xEB),
                intArrayOf(
                    0x44,
                    0x45,
                    0x54,
                    0x55,
                    0x59,
                    0x89,
                    0x95,
                    0x98,
                    0x99,
                    0x9A,
                    0x9D,
                    0xA9,
                    0xAA,
                    0xAE,
                    0xD9,
                    0xEA,
                    0xEE
                ),
                intArrayOf(0x45, 0x55, 0x59, 0x89, 0x95, 0x99, 0x9A, 0x9E, 0xA9, 0xAA, 0xAE, 0xDA, 0xEA, 0xEE),
                intArrayOf(0x45, 0x55, 0x59, 0x5A, 0x95, 0x96, 0x99, 0x9A, 0x9E, 0xA9, 0xAA, 0xAE, 0xDA, 0xEA, 0xEE),
                intArrayOf(
                    0x45,
                    0x55,
                    0x56,
                    0x59,
                    0x5A,
                    0x95,
                    0x96,
                    0x99,
                    0x9A,
                    0x9B,
                    0x9E,
                    0xAA,
                    0xAB,
                    0xAE,
                    0xDA,
                    0xEA,
                    0xEF
                ),
                intArrayOf(0x50, 0x51, 0x54, 0x55, 0x65, 0x91, 0x94, 0x95, 0x96, 0x99, 0xA5, 0xA6, 0xA9, 0xAA, 0xEA),
                intArrayOf(0x51, 0x55, 0x91, 0x95, 0x96, 0x99, 0x9A, 0xA5, 0xA6, 0xA9, 0xAA, 0xE6, 0xEA),
                intArrayOf(0x51, 0x55, 0x56, 0x91, 0x95, 0x96, 0x9A, 0xA5, 0xA6, 0xAA, 0xE6, 0xEA),
                intArrayOf(0x51, 0x55, 0x56, 0x92, 0x95, 0x96, 0x9A, 0xA6, 0xA7, 0xAA, 0xAB, 0xE6, 0xEA, 0xEB),
                intArrayOf(0x54, 0x55, 0x94, 0x95, 0x96, 0x99, 0x9A, 0xA5, 0xA6, 0xA9, 0xAA, 0xE9, 0xEA),
                intArrayOf(0x55, 0x95, 0x96, 0x99, 0x9A, 0xA5, 0xA6, 0xA9, 0xAA, 0xEA),
                intArrayOf(0x55, 0x56, 0x95, 0x96, 0x99, 0x9A, 0xA5, 0xA6, 0xA9, 0xAA, 0xEA),
                intArrayOf(0x55, 0x56, 0x95, 0x96, 0x9A, 0xA6, 0xAA, 0xAB, 0xEA, 0xEB),
                intArrayOf(0x54, 0x55, 0x59, 0x94, 0x95, 0x99, 0x9A, 0xA5, 0xA9, 0xAA, 0xE9, 0xEA),
                intArrayOf(0x55, 0x59, 0x95, 0x96, 0x99, 0x9A, 0xA5, 0xA6, 0xA9, 0xAA, 0xEA),
                intArrayOf(0x45, 0x55, 0x56, 0x59, 0x5A, 0x95, 0x96, 0x99, 0x9A, 0xA5, 0xA6, 0xA9, 0xAA, 0xEA),
                intArrayOf(0x45, 0x55, 0x56, 0x5A, 0x95, 0x96, 0x99, 0x9A, 0xA6, 0xAA, 0xAB, 0xEA, 0xEB),
                intArrayOf(0x54, 0x55, 0x59, 0x95, 0x98, 0x99, 0x9A, 0xA9, 0xAA, 0xAD, 0xAE, 0xE9, 0xEA, 0xEE),
                intArrayOf(0x55, 0x59, 0x95, 0x99, 0x9A, 0xA9, 0xAA, 0xAE, 0xEA, 0xEE),
                intArrayOf(0x45, 0x55, 0x59, 0x5A, 0x95, 0x96, 0x99, 0x9A, 0xA9, 0xAA, 0xAE, 0xEA, 0xEE),
                intArrayOf(0x55, 0x56, 0x59, 0x5A, 0x95, 0x96, 0x99, 0x9A, 0xAA, 0xAB, 0xAE, 0xEA, 0xEF),
                intArrayOf(0x50, 0x51, 0x54, 0x55, 0x65, 0x91, 0x94, 0x95, 0xA5, 0xA6, 0xA9, 0xAA, 0xE5, 0xEA),
                intArrayOf(0x51, 0x55, 0x65, 0x91, 0x95, 0x96, 0xA5, 0xA6, 0xA9, 0xAA, 0xE6, 0xEA),
                intArrayOf(0x51, 0x55, 0x56, 0x65, 0x66, 0x91, 0x95, 0x96, 0xA5, 0xA6, 0xAA, 0xE6, 0xEA),
                intArrayOf(0x51, 0x55, 0x56, 0x66, 0x95, 0x96, 0x9A, 0xA5, 0xA6, 0xA7, 0xAA, 0xAB, 0xE6, 0xEA, 0xEB),
                intArrayOf(0x54, 0x55, 0x65, 0x94, 0x95, 0x99, 0xA5, 0xA6, 0xA9, 0xAA, 0xE9, 0xEA),
                intArrayOf(0x55, 0x65, 0x95, 0x96, 0x99, 0x9A, 0xA5, 0xA6, 0xA9, 0xAA, 0xEA),
                intArrayOf(0x51, 0x55, 0x56, 0x65, 0x66, 0x95, 0x96, 0x99, 0x9A, 0xA5, 0xA6, 0xA9, 0xAA, 0xEA),
                intArrayOf(0x51, 0x55, 0x56, 0x66, 0x95, 0x96, 0x9A, 0xA5, 0xA6, 0xAA, 0xAB, 0xEA, 0xEB),
                intArrayOf(0x54, 0x55, 0x59, 0x65, 0x69, 0x94, 0x95, 0x99, 0xA5, 0xA9, 0xAA, 0xE9, 0xEA),
                intArrayOf(0x54, 0x55, 0x59, 0x65, 0x69, 0x95, 0x96, 0x99, 0x9A, 0xA5, 0xA6, 0xA9, 0xAA, 0xEA),
                intArrayOf(0x55, 0x56, 0x59, 0x65, 0x6A, 0x95, 0x96, 0x99, 0x9A, 0xA5, 0xA6, 0xA9, 0xAA, 0xEA),
                intArrayOf(
                    0x55,
                    0x56,
                    0x5A,
                    0x66,
                    0x6A,
                    0x95,
                    0x96,
                    0x99,
                    0x9A,
                    0xA5,
                    0xA6,
                    0xA9,
                    0xAA,
                    0xAB,
                    0xEA,
                    0xEB
                ),
                intArrayOf(0x54, 0x55, 0x59, 0x69, 0x95, 0x99, 0x9A, 0xA5, 0xA9, 0xAA, 0xAD, 0xAE, 0xE9, 0xEA, 0xEE),
                intArrayOf(0x54, 0x55, 0x59, 0x69, 0x95, 0x99, 0x9A, 0xA5, 0xA9, 0xAA, 0xAE, 0xEA, 0xEE),
                intArrayOf(
                    0x55,
                    0x59,
                    0x5A,
                    0x69,
                    0x6A,
                    0x95,
                    0x96,
                    0x99,
                    0x9A,
                    0xA5,
                    0xA6,
                    0xA9,
                    0xAA,
                    0xAE,
                    0xEA,
                    0xEE
                ),
                intArrayOf(0x55, 0x56, 0x59, 0x5A, 0x6A, 0x95, 0x96, 0x99, 0x9A, 0xA6, 0xA9, 0xAA, 0xAB, 0xAE, 0xEA),
                intArrayOf(
                    0x50,
                    0x51,
                    0x54,
                    0x55,
                    0x65,
                    0x95,
                    0xA1,
                    0xA4,
                    0xA5,
                    0xA6,
                    0xA9,
                    0xAA,
                    0xB5,
                    0xBA,
                    0xE5,
                    0xEA,
                    0xFA
                ),
                intArrayOf(0x51, 0x55, 0x65, 0x95, 0xA1, 0xA5, 0xA6, 0xA9, 0xAA, 0xB6, 0xBA, 0xE6, 0xEA, 0xFA),
                intArrayOf(0x51, 0x55, 0x65, 0x66, 0x95, 0x96, 0xA5, 0xA6, 0xA9, 0xAA, 0xB6, 0xBA, 0xE6, 0xEA, 0xFA),
                intArrayOf(
                    0x51,
                    0x55,
                    0x56,
                    0x65,
                    0x66,
                    0x95,
                    0x96,
                    0xA5,
                    0xA6,
                    0xA7,
                    0xAA,
                    0xAB,
                    0xB6,
                    0xBA,
                    0xE6,
                    0xEA,
                    0xFB
                ),
                intArrayOf(0x54, 0x55, 0x65, 0x95, 0xA4, 0xA5, 0xA6, 0xA9, 0xAA, 0xB9, 0xBA, 0xE9, 0xEA, 0xFA),
                intArrayOf(0x55, 0x65, 0x95, 0xA5, 0xA6, 0xA9, 0xAA, 0xBA, 0xEA, 0xFA),
                intArrayOf(0x51, 0x55, 0x65, 0x66, 0x95, 0x96, 0xA5, 0xA6, 0xA9, 0xAA, 0xBA, 0xEA, 0xFA),
                intArrayOf(0x55, 0x56, 0x65, 0x66, 0x95, 0x96, 0xA5, 0xA6, 0xAA, 0xAB, 0xBA, 0xEA, 0xFB),
                intArrayOf(0x54, 0x55, 0x65, 0x69, 0x95, 0x99, 0xA5, 0xA6, 0xA9, 0xAA, 0xB9, 0xBA, 0xE9, 0xEA, 0xFA),
                intArrayOf(0x54, 0x55, 0x65, 0x69, 0x95, 0x99, 0xA5, 0xA6, 0xA9, 0xAA, 0xBA, 0xEA, 0xFA),
                intArrayOf(
                    0x55,
                    0x65,
                    0x66,
                    0x69,
                    0x6A,
                    0x95,
                    0x96,
                    0x99,
                    0x9A,
                    0xA5,
                    0xA6,
                    0xA9,
                    0xAA,
                    0xBA,
                    0xEA,
                    0xFA
                ),
                intArrayOf(0x55, 0x56, 0x65, 0x66, 0x6A, 0x95, 0x96, 0x9A, 0xA5, 0xA6, 0xA9, 0xAA, 0xAB, 0xBA, 0xEA),
                intArrayOf(
                    0x54,
                    0x55,
                    0x59,
                    0x65,
                    0x69,
                    0x95,
                    0x99,
                    0xA5,
                    0xA9,
                    0xAA,
                    0xAD,
                    0xAE,
                    0xB9,
                    0xBA,
                    0xE9,
                    0xEA,
                    0xFE
                ),
                intArrayOf(0x55, 0x59, 0x65, 0x69, 0x95, 0x99, 0xA5, 0xA9, 0xAA, 0xAE, 0xBA, 0xEA, 0xFE),
                intArrayOf(0x55, 0x59, 0x65, 0x69, 0x6A, 0x95, 0x99, 0x9A, 0xA5, 0xA6, 0xA9, 0xAA, 0xAE, 0xBA, 0xEA),
                intArrayOf(
                    0x55,
                    0x56,
                    0x59,
                    0x5A,
                    0x65,
                    0x66,
                    0x69,
                    0x6A,
                    0x95,
                    0x96,
                    0x99,
                    0x9A,
                    0xA5,
                    0xA6,
                    0xA9,
                    0xAA,
                    0xAB,
                    0xAE,
                    0xBA,
                    0xEA
                )
            )
            val latticePoints = arrayOfNulls<LatticePoint4D>(256)
            for (i in 0..255) {
                val cx = (i shr 0 and 3) - 1
                val cy = (i shr 2 and 3) - 1
                val cz = (i shr 4 and 3) - 1
                val cw = (i shr 6 and 3) - 1
                latticePoints[i] = LatticePoint4D(cx, cy, cz, cw)
            }
            for (i in 0..255) {
                LOOKUP_4D[i] = arrayOfNulls(lookup4DPregen[i].size)
                for (j in lookup4DPregen[i].indices) {
                    LOOKUP_4D[i]?.set(j, latticePoints[lookup4DPregen[i][j]])
                }
            }
        }

        init {
            GRADIENTS_2D = arrayOfNulls(PSIZE)
            val grad2 = arrayOf(
                Grad2(0.130526192220052, 0.99144486137381),
                Grad2(0.38268343236509, 0.923879532511287),
                Grad2(0.608761429008721, 0.793353340291235),
                Grad2(0.793353340291235, 0.608761429008721),
                Grad2(0.923879532511287, 0.38268343236509),
                Grad2(0.99144486137381, 0.130526192220051),
                Grad2(0.99144486137381, -0.130526192220051),
                Grad2(0.923879532511287, -0.38268343236509),
                Grad2(0.793353340291235, -0.60876142900872),
                Grad2(0.608761429008721, -0.793353340291235),
                Grad2(0.38268343236509, -0.923879532511287),
                Grad2(0.130526192220052, -0.99144486137381),
                Grad2(-0.130526192220052, -0.99144486137381),
                Grad2(-0.38268343236509, -0.923879532511287),
                Grad2(-0.608761429008721, -0.793353340291235),
                Grad2(-0.793353340291235, -0.608761429008721),
                Grad2(-0.923879532511287, -0.38268343236509),
                Grad2(-0.99144486137381, -0.130526192220052),
                Grad2(-0.99144486137381, 0.130526192220051),
                Grad2(-0.923879532511287, 0.38268343236509),
                Grad2(-0.793353340291235, 0.608761429008721),
                Grad2(-0.608761429008721, 0.793353340291235),
                Grad2(-0.38268343236509, 0.923879532511287),
                Grad2(-0.130526192220052, 0.99144486137381)
            )
            val grad2XBeforeY = arrayOfNulls<Grad2>(grad2.size)
            for (i in grad2.indices) {
                grad2[i].dx /= N2
                grad2[i].dy /= N2
            }
            for (i in 0 until PSIZE) {
                GRADIENTS_2D[i] = grad2[i % grad2.size]
            }
            GRADIENTS_3D = arrayOfNulls(PSIZE)
            val grad3 = arrayOf(
                Grad3(-2.22474487139, -2.22474487139, -1.0),
                Grad3(-2.22474487139, -2.22474487139, 1.0),
                Grad3(-3.0862664687972017, -1.1721513422464978, 0.0),
                Grad3(-1.1721513422464978, -3.0862664687972017, 0.0),
                Grad3(-2.22474487139, -1.0, -2.22474487139),
                Grad3(-2.22474487139, 1.0, -2.22474487139),
                Grad3(-1.1721513422464978, 0.0, -3.0862664687972017),
                Grad3(-3.0862664687972017, 0.0, -1.1721513422464978),
                Grad3(-2.22474487139, -1.0, 2.22474487139),
                Grad3(-2.22474487139, 1.0, 2.22474487139),
                Grad3(-3.0862664687972017, 0.0, 1.1721513422464978),
                Grad3(-1.1721513422464978, 0.0, 3.0862664687972017),
                Grad3(-2.22474487139, 2.22474487139, -1.0),
                Grad3(-2.22474487139, 2.22474487139, 1.0),
                Grad3(-1.1721513422464978, 3.0862664687972017, 0.0),
                Grad3(-3.0862664687972017, 1.1721513422464978, 0.0),
                Grad3(-1.0, -2.22474487139, -2.22474487139),
                Grad3(1.0, -2.22474487139, -2.22474487139),
                Grad3(0.0, -3.0862664687972017, -1.1721513422464978),
                Grad3(0.0, -1.1721513422464978, -3.0862664687972017),
                Grad3(-1.0, -2.22474487139, 2.22474487139),
                Grad3(1.0, -2.22474487139, 2.22474487139),
                Grad3(0.0, -1.1721513422464978, 3.0862664687972017),
                Grad3(0.0, -3.0862664687972017, 1.1721513422464978),
                Grad3(-1.0, 2.22474487139, -2.22474487139),
                Grad3(1.0, 2.22474487139, -2.22474487139),
                Grad3(0.0, 1.1721513422464978, -3.0862664687972017),
                Grad3(0.0, 3.0862664687972017, -1.1721513422464978),
                Grad3(-1.0, 2.22474487139, 2.22474487139),
                Grad3(1.0, 2.22474487139, 2.22474487139),
                Grad3(0.0, 3.0862664687972017, 1.1721513422464978),
                Grad3(0.0, 1.1721513422464978, 3.0862664687972017),
                Grad3(2.22474487139, -2.22474487139, -1.0),
                Grad3(2.22474487139, -2.22474487139, 1.0),
                Grad3(1.1721513422464978, -3.0862664687972017, 0.0),
                Grad3(3.0862664687972017, -1.1721513422464978, 0.0),
                Grad3(2.22474487139, -1.0, -2.22474487139),
                Grad3(2.22474487139, 1.0, -2.22474487139),
                Grad3(3.0862664687972017, 0.0, -1.1721513422464978),
                Grad3(1.1721513422464978, 0.0, -3.0862664687972017),
                Grad3(2.22474487139, -1.0, 2.22474487139),
                Grad3(2.22474487139, 1.0, 2.22474487139),
                Grad3(1.1721513422464978, 0.0, 3.0862664687972017),
                Grad3(3.0862664687972017, 0.0, 1.1721513422464978),
                Grad3(2.22474487139, 2.22474487139, -1.0),
                Grad3(2.22474487139, 2.22474487139, 1.0),
                Grad3(3.0862664687972017, 1.1721513422464978, 0.0),
                Grad3(1.1721513422464978, 3.0862664687972017, 0.0)
            )
            for (i in grad3.indices) {
                grad3[i].dx /= N3
                grad3[i].dy /= N3
                grad3[i].dz /= N3
            }
            for (i in 0 until PSIZE) {
                GRADIENTS_3D[i] = grad3[i % grad3.size]
            }
            GRADIENTS_4D = arrayOfNulls(PSIZE)
            val grad4 = arrayOf(
                Grad4(-0.753341017856078, -0.37968289875261624, -0.37968289875261624, -0.37968289875261624),
                Grad4(-0.7821684431180708, -0.4321472685365301, -0.4321472685365301, 0.12128480194602098),
                Grad4(-0.7821684431180708, -0.4321472685365301, 0.12128480194602098, -0.4321472685365301),
                Grad4(-0.7821684431180708, 0.12128480194602098, -0.4321472685365301, -0.4321472685365301),
                Grad4(-0.8586508742123365, -0.508629699630796, 0.044802370851755174, 0.044802370851755174),
                Grad4(-0.8586508742123365, 0.044802370851755174, -0.508629699630796, 0.044802370851755174),
                Grad4(-0.8586508742123365, 0.044802370851755174, 0.044802370851755174, -0.508629699630796),
                Grad4(-0.9982828964265062, -0.03381941603233842, -0.03381941603233842, -0.03381941603233842),
                Grad4(-0.37968289875261624, -0.753341017856078, -0.37968289875261624, -0.37968289875261624),
                Grad4(-0.4321472685365301, -0.7821684431180708, -0.4321472685365301, 0.12128480194602098),
                Grad4(-0.4321472685365301, -0.7821684431180708, 0.12128480194602098, -0.4321472685365301),
                Grad4(0.12128480194602098, -0.7821684431180708, -0.4321472685365301, -0.4321472685365301),
                Grad4(-0.508629699630796, -0.8586508742123365, 0.044802370851755174, 0.044802370851755174),
                Grad4(0.044802370851755174, -0.8586508742123365, -0.508629699630796, 0.044802370851755174),
                Grad4(0.044802370851755174, -0.8586508742123365, 0.044802370851755174, -0.508629699630796),
                Grad4(-0.03381941603233842, -0.9982828964265062, -0.03381941603233842, -0.03381941603233842),
                Grad4(-0.37968289875261624, -0.37968289875261624, -0.753341017856078, -0.37968289875261624),
                Grad4(-0.4321472685365301, -0.4321472685365301, -0.7821684431180708, 0.12128480194602098),
                Grad4(-0.4321472685365301, 0.12128480194602098, -0.7821684431180708, -0.4321472685365301),
                Grad4(0.12128480194602098, -0.4321472685365301, -0.7821684431180708, -0.4321472685365301),
                Grad4(-0.508629699630796, 0.044802370851755174, -0.8586508742123365, 0.044802370851755174),
                Grad4(0.044802370851755174, -0.508629699630796, -0.8586508742123365, 0.044802370851755174),
                Grad4(0.044802370851755174, 0.044802370851755174, -0.8586508742123365, -0.508629699630796),
                Grad4(-0.03381941603233842, -0.03381941603233842, -0.9982828964265062, -0.03381941603233842),
                Grad4(-0.37968289875261624, -0.37968289875261624, -0.37968289875261624, -0.753341017856078),
                Grad4(-0.4321472685365301, -0.4321472685365301, 0.12128480194602098, -0.7821684431180708),
                Grad4(-0.4321472685365301, 0.12128480194602098, -0.4321472685365301, -0.7821684431180708),
                Grad4(0.12128480194602098, -0.4321472685365301, -0.4321472685365301, -0.7821684431180708),
                Grad4(-0.508629699630796, 0.044802370851755174, 0.044802370851755174, -0.8586508742123365),
                Grad4(0.044802370851755174, -0.508629699630796, 0.044802370851755174, -0.8586508742123365),
                Grad4(0.044802370851755174, 0.044802370851755174, -0.508629699630796, -0.8586508742123365),
                Grad4(-0.03381941603233842, -0.03381941603233842, -0.03381941603233842, -0.9982828964265062),
                Grad4(-0.6740059517812944, -0.3239847771997537, -0.3239847771997537, 0.5794684678643381),
                Grad4(-0.7504883828755602, -0.4004672082940195, 0.15296486218853164, 0.5029860367700724),
                Grad4(-0.7504883828755602, 0.15296486218853164, -0.4004672082940195, 0.5029860367700724),
                Grad4(-0.8828161875373585, 0.08164729285680945, 0.08164729285680945, 0.4553054119602712),
                Grad4(-0.4553054119602712, -0.08164729285680945, -0.08164729285680945, 0.8828161875373585),
                Grad4(-0.5029860367700724, -0.15296486218853164, 0.4004672082940195, 0.7504883828755602),
                Grad4(-0.5029860367700724, 0.4004672082940195, -0.15296486218853164, 0.7504883828755602),
                Grad4(-0.5794684678643381, 0.3239847771997537, 0.3239847771997537, 0.6740059517812944),
                Grad4(-0.3239847771997537, -0.6740059517812944, -0.3239847771997537, 0.5794684678643381),
                Grad4(-0.4004672082940195, -0.7504883828755602, 0.15296486218853164, 0.5029860367700724),
                Grad4(0.15296486218853164, -0.7504883828755602, -0.4004672082940195, 0.5029860367700724),
                Grad4(0.08164729285680945, -0.8828161875373585, 0.08164729285680945, 0.4553054119602712),
                Grad4(-0.08164729285680945, -0.4553054119602712, -0.08164729285680945, 0.8828161875373585),
                Grad4(-0.15296486218853164, -0.5029860367700724, 0.4004672082940195, 0.7504883828755602),
                Grad4(0.4004672082940195, -0.5029860367700724, -0.15296486218853164, 0.7504883828755602),
                Grad4(0.3239847771997537, -0.5794684678643381, 0.3239847771997537, 0.6740059517812944),
                Grad4(-0.3239847771997537, -0.3239847771997537, -0.6740059517812944, 0.5794684678643381),
                Grad4(-0.4004672082940195, 0.15296486218853164, -0.7504883828755602, 0.5029860367700724),
                Grad4(0.15296486218853164, -0.4004672082940195, -0.7504883828755602, 0.5029860367700724),
                Grad4(0.08164729285680945, 0.08164729285680945, -0.8828161875373585, 0.4553054119602712),
                Grad4(-0.08164729285680945, -0.08164729285680945, -0.4553054119602712, 0.8828161875373585),
                Grad4(-0.15296486218853164, 0.4004672082940195, -0.5029860367700724, 0.7504883828755602),
                Grad4(0.4004672082940195, -0.15296486218853164, -0.5029860367700724, 0.7504883828755602),
                Grad4(0.3239847771997537, 0.3239847771997537, -0.5794684678643381, 0.6740059517812944),
                Grad4(-0.6740059517812944, -0.3239847771997537, 0.5794684678643381, -0.3239847771997537),
                Grad4(-0.7504883828755602, -0.4004672082940195, 0.5029860367700724, 0.15296486218853164),
                Grad4(-0.7504883828755602, 0.15296486218853164, 0.5029860367700724, -0.4004672082940195),
                Grad4(-0.8828161875373585, 0.08164729285680945, 0.4553054119602712, 0.08164729285680945),
                Grad4(-0.4553054119602712, -0.08164729285680945, 0.8828161875373585, -0.08164729285680945),
                Grad4(-0.5029860367700724, -0.15296486218853164, 0.7504883828755602, 0.4004672082940195),
                Grad4(-0.5029860367700724, 0.4004672082940195, 0.7504883828755602, -0.15296486218853164),
                Grad4(-0.5794684678643381, 0.3239847771997537, 0.6740059517812944, 0.3239847771997537),
                Grad4(-0.3239847771997537, -0.6740059517812944, 0.5794684678643381, -0.3239847771997537),
                Grad4(-0.4004672082940195, -0.7504883828755602, 0.5029860367700724, 0.15296486218853164),
                Grad4(0.15296486218853164, -0.7504883828755602, 0.5029860367700724, -0.4004672082940195),
                Grad4(0.08164729285680945, -0.8828161875373585, 0.4553054119602712, 0.08164729285680945),
                Grad4(-0.08164729285680945, -0.4553054119602712, 0.8828161875373585, -0.08164729285680945),
                Grad4(-0.15296486218853164, -0.5029860367700724, 0.7504883828755602, 0.4004672082940195),
                Grad4(0.4004672082940195, -0.5029860367700724, 0.7504883828755602, -0.15296486218853164),
                Grad4(0.3239847771997537, -0.5794684678643381, 0.6740059517812944, 0.3239847771997537),
                Grad4(-0.3239847771997537, -0.3239847771997537, 0.5794684678643381, -0.6740059517812944),
                Grad4(-0.4004672082940195, 0.15296486218853164, 0.5029860367700724, -0.7504883828755602),
                Grad4(0.15296486218853164, -0.4004672082940195, 0.5029860367700724, -0.7504883828755602),
                Grad4(0.08164729285680945, 0.08164729285680945, 0.4553054119602712, -0.8828161875373585),
                Grad4(-0.08164729285680945, -0.08164729285680945, 0.8828161875373585, -0.4553054119602712),
                Grad4(-0.15296486218853164, 0.4004672082940195, 0.7504883828755602, -0.5029860367700724),
                Grad4(0.4004672082940195, -0.15296486218853164, 0.7504883828755602, -0.5029860367700724),
                Grad4(0.3239847771997537, 0.3239847771997537, 0.6740059517812944, -0.5794684678643381),
                Grad4(-0.6740059517812944, 0.5794684678643381, -0.3239847771997537, -0.3239847771997537),
                Grad4(-0.7504883828755602, 0.5029860367700724, -0.4004672082940195, 0.15296486218853164),
                Grad4(-0.7504883828755602, 0.5029860367700724, 0.15296486218853164, -0.4004672082940195),
                Grad4(-0.8828161875373585, 0.4553054119602712, 0.08164729285680945, 0.08164729285680945),
                Grad4(-0.4553054119602712, 0.8828161875373585, -0.08164729285680945, -0.08164729285680945),
                Grad4(-0.5029860367700724, 0.7504883828755602, -0.15296486218853164, 0.4004672082940195),
                Grad4(-0.5029860367700724, 0.7504883828755602, 0.4004672082940195, -0.15296486218853164),
                Grad4(-0.5794684678643381, 0.6740059517812944, 0.3239847771997537, 0.3239847771997537),
                Grad4(-0.3239847771997537, 0.5794684678643381, -0.6740059517812944, -0.3239847771997537),
                Grad4(-0.4004672082940195, 0.5029860367700724, -0.7504883828755602, 0.15296486218853164),
                Grad4(0.15296486218853164, 0.5029860367700724, -0.7504883828755602, -0.4004672082940195),
                Grad4(0.08164729285680945, 0.4553054119602712, -0.8828161875373585, 0.08164729285680945),
                Grad4(-0.08164729285680945, 0.8828161875373585, -0.4553054119602712, -0.08164729285680945),
                Grad4(-0.15296486218853164, 0.7504883828755602, -0.5029860367700724, 0.4004672082940195),
                Grad4(0.4004672082940195, 0.7504883828755602, -0.5029860367700724, -0.15296486218853164),
                Grad4(0.3239847771997537, 0.6740059517812944, -0.5794684678643381, 0.3239847771997537),
                Grad4(-0.3239847771997537, 0.5794684678643381, -0.3239847771997537, -0.6740059517812944),
                Grad4(-0.4004672082940195, 0.5029860367700724, 0.15296486218853164, -0.7504883828755602),
                Grad4(0.15296486218853164, 0.5029860367700724, -0.4004672082940195, -0.7504883828755602),
                Grad4(0.08164729285680945, 0.4553054119602712, 0.08164729285680945, -0.8828161875373585),
                Grad4(-0.08164729285680945, 0.8828161875373585, -0.08164729285680945, -0.4553054119602712),
                Grad4(-0.15296486218853164, 0.7504883828755602, 0.4004672082940195, -0.5029860367700724),
                Grad4(0.4004672082940195, 0.7504883828755602, -0.15296486218853164, -0.5029860367700724),
                Grad4(0.3239847771997537, 0.6740059517812944, 0.3239847771997537, -0.5794684678643381),
                Grad4(0.5794684678643381, -0.6740059517812944, -0.3239847771997537, -0.3239847771997537),
                Grad4(0.5029860367700724, -0.7504883828755602, -0.4004672082940195, 0.15296486218853164),
                Grad4(0.5029860367700724, -0.7504883828755602, 0.15296486218853164, -0.4004672082940195),
                Grad4(0.4553054119602712, -0.8828161875373585, 0.08164729285680945, 0.08164729285680945),
                Grad4(0.8828161875373585, -0.4553054119602712, -0.08164729285680945, -0.08164729285680945),
                Grad4(0.7504883828755602, -0.5029860367700724, -0.15296486218853164, 0.4004672082940195),
                Grad4(0.7504883828755602, -0.5029860367700724, 0.4004672082940195, -0.15296486218853164),
                Grad4(0.6740059517812944, -0.5794684678643381, 0.3239847771997537, 0.3239847771997537),
                Grad4(0.5794684678643381, -0.3239847771997537, -0.6740059517812944, -0.3239847771997537),
                Grad4(0.5029860367700724, -0.4004672082940195, -0.7504883828755602, 0.15296486218853164),
                Grad4(0.5029860367700724, 0.15296486218853164, -0.7504883828755602, -0.4004672082940195),
                Grad4(0.4553054119602712, 0.08164729285680945, -0.8828161875373585, 0.08164729285680945),
                Grad4(0.8828161875373585, -0.08164729285680945, -0.4553054119602712, -0.08164729285680945),
                Grad4(0.7504883828755602, -0.15296486218853164, -0.5029860367700724, 0.4004672082940195),
                Grad4(0.7504883828755602, 0.4004672082940195, -0.5029860367700724, -0.15296486218853164),
                Grad4(0.6740059517812944, 0.3239847771997537, -0.5794684678643381, 0.3239847771997537),
                Grad4(0.5794684678643381, -0.3239847771997537, -0.3239847771997537, -0.6740059517812944),
                Grad4(0.5029860367700724, -0.4004672082940195, 0.15296486218853164, -0.7504883828755602),
                Grad4(0.5029860367700724, 0.15296486218853164, -0.4004672082940195, -0.7504883828755602),
                Grad4(0.4553054119602712, 0.08164729285680945, 0.08164729285680945, -0.8828161875373585),
                Grad4(0.8828161875373585, -0.08164729285680945, -0.08164729285680945, -0.4553054119602712),
                Grad4(0.7504883828755602, -0.15296486218853164, 0.4004672082940195, -0.5029860367700724),
                Grad4(0.7504883828755602, 0.4004672082940195, -0.15296486218853164, -0.5029860367700724),
                Grad4(0.6740059517812944, 0.3239847771997537, 0.3239847771997537, -0.5794684678643381),
                Grad4(0.03381941603233842, 0.03381941603233842, 0.03381941603233842, 0.9982828964265062),
                Grad4(-0.044802370851755174, -0.044802370851755174, 0.508629699630796, 0.8586508742123365),
                Grad4(-0.044802370851755174, 0.508629699630796, -0.044802370851755174, 0.8586508742123365),
                Grad4(-0.12128480194602098, 0.4321472685365301, 0.4321472685365301, 0.7821684431180708),
                Grad4(0.508629699630796, -0.044802370851755174, -0.044802370851755174, 0.8586508742123365),
                Grad4(0.4321472685365301, -0.12128480194602098, 0.4321472685365301, 0.7821684431180708),
                Grad4(0.4321472685365301, 0.4321472685365301, -0.12128480194602098, 0.7821684431180708),
                Grad4(0.37968289875261624, 0.37968289875261624, 0.37968289875261624, 0.753341017856078),
                Grad4(0.03381941603233842, 0.03381941603233842, 0.9982828964265062, 0.03381941603233842),
                Grad4(-0.044802370851755174, 0.044802370851755174, 0.8586508742123365, 0.508629699630796),
                Grad4(-0.044802370851755174, 0.508629699630796, 0.8586508742123365, -0.044802370851755174),
                Grad4(-0.12128480194602098, 0.4321472685365301, 0.7821684431180708, 0.4321472685365301),
                Grad4(0.508629699630796, -0.044802370851755174, 0.8586508742123365, -0.044802370851755174),
                Grad4(0.4321472685365301, -0.12128480194602098, 0.7821684431180708, 0.4321472685365301),
                Grad4(0.4321472685365301, 0.4321472685365301, 0.7821684431180708, -0.12128480194602098),
                Grad4(0.37968289875261624, 0.37968289875261624, 0.753341017856078, 0.37968289875261624),
                Grad4(0.03381941603233842, 0.9982828964265062, 0.03381941603233842, 0.03381941603233842),
                Grad4(-0.044802370851755174, 0.8586508742123365, -0.044802370851755174, 0.508629699630796),
                Grad4(-0.044802370851755174, 0.8586508742123365, 0.508629699630796, -0.044802370851755174),
                Grad4(-0.12128480194602098, 0.7821684431180708, 0.4321472685365301, 0.4321472685365301),
                Grad4(0.508629699630796, 0.8586508742123365, -0.044802370851755174, -0.044802370851755174),
                Grad4(0.4321472685365301, 0.7821684431180708, -0.12128480194602098, 0.4321472685365301),
                Grad4(0.4321472685365301, 0.7821684431180708, 0.4321472685365301, -0.12128480194602098),
                Grad4(0.37968289875261624, 0.753341017856078, 0.37968289875261624, 0.37968289875261624),
                Grad4(0.9982828964265062, 0.03381941603233842, 0.03381941603233842, 0.03381941603233842),
                Grad4(0.8586508742123365, -0.044802370851755174, -0.044802370851755174, 0.508629699630796),
                Grad4(0.8586508742123365, -0.044802370851755174, 0.508629699630796, -0.044802370851755174),
                Grad4(0.7821684431180708, -0.12128480194602098, 0.4321472685365301, 0.4321472685365301),
                Grad4(0.8586508742123365, 0.508629699630796, -0.044802370851755174, -0.044802370851755174),
                Grad4(0.7821684431180708, 0.4321472685365301, -0.12128480194602098, 0.4321472685365301),
                Grad4(0.7821684431180708, 0.4321472685365301, 0.4321472685365301, -0.12128480194602098),
                Grad4(0.753341017856078, 0.37968289875261624, 0.37968289875261624, 0.37968289875261624)
            )
            for (i in grad4.indices) {
                grad4[i].dx /= N4
                grad4[i].dy /= N4
                grad4[i].dz /= N4
                grad4[i].dw /= N4
            }
            for (i in 0 until PSIZE) {
                GRADIENTS_4D[i] = grad4[i % grad4.size]
            }
        }
    }

    private class LatticePoint2D(var xsv: Int, var ysv: Int) {
        var dx: Double
        var dy: Double

        init {
            val ssv = (xsv + ysv) * -0.211324865405187
            dx = -xsv - ssv
            dy = -ysv - ssv
        }
    }

    private class LatticePoint3D(xrv: Int, yrv: Int, zrv: Int, lattice: Int) {
        var dxr: Double
        var dyr: Double
        var dzr: Double
        var xrv: Int
        var yrv: Int
        var zrv: Int
        var nextOnFailure: LatticePoint3D? = null
        var nextOnSuccess: LatticePoint3D? = null

        init {
            dxr = -xrv + lattice * 0.5
            dyr = -yrv + lattice * 0.5
            dzr = -zrv + lattice * 0.5
            this.xrv = xrv + lattice * 1024
            this.yrv = yrv + lattice * 1024
            this.zrv = zrv + lattice * 1024
        }
    }

    private class LatticePoint4D(var xsv: Int, var ysv: Int, var zsv: Int, var wsv: Int) {
        var dx: Double
        var dy: Double
        var dz: Double
        var dw: Double

        init {
            val ssv = (xsv + ysv + zsv + wsv) * -0.138196601125011
            dx = -xsv - ssv
            dy = -ysv - ssv
            dz = -zsv - ssv
            dw = -wsv - ssv
        }
    }

    /*
	 * Gradients
	 */
    class Grad2(var dx: Double, var dy: Double)

    class Grad3(var dx: Double, var dy: Double, var dz: Double)

    class Grad4(var dx: Double, var dy: Double, var dz: Double, var dw: Double)

    init {
        var seed = seed
        perm = ShortArray(PSIZE)
        permGrad2 = arrayOfNulls(PSIZE)
        permGrad3 = arrayOfNulls(PSIZE)
        permGrad4 = arrayOfNulls(PSIZE)
        val source = ShortArray(PSIZE)
        for (i in 0 until PSIZE) source[i] = i.toShort()
        for (i in PSIZE - 1 downTo 0) {
            seed = seed * 6364136223846793005L + 1442695040888963407L
            var r = ((seed + 31) % (i + 1)).toInt()
            if (r < 0) r += i + 1
            perm[i] = source[r]
            permGrad2[i] = GRADIENTS_2D[perm[i].toInt()]
            permGrad3[i] = GRADIENTS_3D[perm[i].toInt()]
            permGrad4[i] = GRADIENTS_4D[perm[i].toInt()]
            source[r] = source[i]
        }
    }
}