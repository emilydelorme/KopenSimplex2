package opensimplex

import kotlin.experimental.xor

/**
 * Converted to kotlin from https://github.com/KdotJPG/OpenSimplex2
 *
 * K.jpg's OpenSimplex 2, faster variant
 *
 * - 2D is standard simplex implemented using a lookup table.
 * - 3D is "Re-oriented 4-point BCC noise" which constructs an
 * isomorphic BCC lattice in a much different way than usual.
 *
 * Multiple versions of each function are provided. See the
 * documentation above each, for more info.
 */
class OpenSimplex2F(seed: Long) {
    private val perm: ShortArray
    private val permGrad2: Array<Grad2?>
    private val permGrad3: Array<Grad3?>
    /*
	 * Noise Evaluators
	 */
    /**
     * 2D Simplex noise, standard lattice orientation.
     */
    fun noise2(x: Double, y: Double): Double {

        // Get points for A2* lattice
        val s = 0.366025403784439 * (x + y)
        val xs = x + s
        val ys = y + s
        return noise2_Base(xs, ys)
    }

    /**
     * 2D Simplex noise, with Y pointing down the main diagonal.
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
     * 2D Simplex noise base.
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
        val index = ((ysi - xsi) / 2 + 1).toInt() * 3
        val ssi = (xsi + ysi) * -0.211324865405187
        val xi = xsi + ssi
        val yi = ysi + ssi

        // Point contributions
        for (i in 0..2) {
            val c = LOOKUP_2D[index + i]
            val dx = xi + c!!.dx
            val dy = yi + c.dy
            var attn = 0.5 - dx * dx - dy * dy
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
     * 3D Re-oriented 4-point BCC noise, classic orientation.
     * Proper substitute for 3D Simplex in light of Forbidden Formulae.
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
     * 3D Re-oriented 4-point BCC noise, with better visual isotropy in (X, Y).
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
     * 3D Re-oriented 4-point BCC noise, with better visual isotropy in (X, Z).
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
     * than to build up the index with enough info to isolate 4 points.
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
            var attn = 0.5 - dxr * dxr - dyr * dyr - dzr * dzr
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
        const val N2 = 0.01001634121365712
        const val N3 = 0.030485933181293584
        private val GRADIENTS_2D: Array<Grad2?>
        private val GRADIENTS_3D: Array<Grad3?>

        init {
            LOOKUP_2D = arrayOfNulls(2 * 3)
            LOOKUP_3D = arrayOfNulls(8)
            for (i in 0..1) {
                var i1: Int
                var j1: Int
                if (i and 1 == 0) {
                    i1 = 1
                    j1 = 0
                } else {
                    i1 = 0
                    j1 = 1
                }
                LOOKUP_2D[i * 3 + 0] = LatticePoint2D(0, 0)
                LOOKUP_2D[i * 3 + 1] = LatticePoint2D(1, 1)
                LOOKUP_2D[i * 3 + 2] = LatticePoint2D(i1, j1)
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

                // Each single step away on the first half-lattice.
                val c2 = LatticePoint3D(i1 xor 1, j1, k1, 0)
                val c3 = LatticePoint3D(i1, j1 xor 1, k1, 0)
                val c4 = LatticePoint3D(i1, j1, k1 xor 1, 0)

                // Each single step away on the second half-lattice.
                val c5 =
                    LatticePoint3D(i1 + (i2 xor 1), j1 + j2, k1 + k2, 1)
                val c6 =
                    LatticePoint3D(i1 + i2, j1 + (j2 xor 1), k1 + k2, 1)
                val c7 =
                    LatticePoint3D(i1 + i2, j1 + j2, k1 + (k2 xor 1), 1)

                // First two are guaranteed.
                c0.nextOnSuccess = c1
                c0.nextOnFailure = c0.nextOnSuccess
                c1.nextOnSuccess = c2
                c1.nextOnFailure = c1.nextOnSuccess

                // Once we find one on the first half-lattice, the rest are out.
                // In addition, knowing c2 rules out c5.
                c2.nextOnFailure = c3
                c2.nextOnSuccess = c6
                c3.nextOnFailure = c4
                c3.nextOnSuccess = c5
                c4.nextOnSuccess = c5
                c4.nextOnFailure = c4.nextOnSuccess

                // Once we find one on the second half-lattice, the rest are out.
                c5.nextOnFailure = c6
                c5.nextOnSuccess = null
                c6.nextOnFailure = c7
                c6.nextOnSuccess = null
                c7.nextOnSuccess = null
                c7.nextOnFailure = c7.nextOnSuccess
                LOOKUP_3D[i] = c0
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

    /*
	 * Gradients
	 */
    class Grad2(var dx: Double, var dy: Double)

    class Grad3(var dx: Double, var dy: Double, var dz: Double)

    init {
        var seed = seed
        perm = ShortArray(PSIZE)
        permGrad2 = arrayOfNulls(PSIZE)
        permGrad3 = arrayOfNulls(PSIZE)
        val source = ShortArray(PSIZE)
        for (i in 0 until PSIZE) source[i] = i.toShort()
        for (i in PSIZE - 1 downTo 0) {
            seed = seed * 6364136223846793005L + 1442695040888963407L
            var r = ((seed + 31) % (i + 1)).toInt()
            if (r < 0) r += i + 1
            perm[i] = source[r]
            permGrad2[i] = GRADIENTS_2D[perm[i].toInt()]
            permGrad3[i] = GRADIENTS_3D[perm[i].toInt()]
            source[r] = source[i]
        }
    }
}
