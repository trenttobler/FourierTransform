using System;
using System.Linq;
using System.Numerics;

namespace TrentTobler.Algorithms.FourierTransform
{
	public static class FastFourierTransform
	{
		/// <summary>
		/// Perform a fast fourier transform.
		/// </summary>
		/// <param name="data">The data to transform.  It must be an exact power of 2 in length.</param>
		public static void Fft( this Complex[] data )
		{
			var n = data.Length;
			var log2N = FloorLog2( n );
			if( ( 1 << log2N ) != n )
				throw new NotSupportedException( "only support FFT on array length of exact power of 2" );

			FftMultiDPlane( data, 0, 0, log2N );
		}

		/// <summary>
		/// Perform an inverse fast fourier transform.
		/// </summary>
		/// <param name="data">The data to transform.  It must be an exact power of 2 in length.</param>
		public static void InverseFft( this Complex[] data )
		{
			var n = data.Length;
			var log2N = FloorLog2( n );
			if( ( 1 << log2N ) != n )
				throw new NotSupportedException( "only support FFT on array length of exact power of 2" );

			Array.Reverse( data, 1, n - 1 );
			data.Fft();

			for( var i = 0; i < n; ++i )
			{
				data[i] /= n;
			}
		}

		/// <summary>
		/// Perform a multi-dimensional fast fourier transform.
		/// </summary>
		/// <param name="data">The data to transform.</param>
		/// <param name="log2ns">The log2 order on each dimension.</param>
		public static void MultiFft( this Complex[] data, params byte[] log2ns )
		{
			var dn = log2ns.Length;
			var totalBits = log2ns.Select( b => (int) b ).Sum();

			if( data.Length < ( 1 << totalBits ) )
			{
				throw new ArgumentException( "dimensional mismatch" );
			}

			var lowerBits = 0;
			var lowerMask = 0;
			for( var d = 0; d < dn; ++d )
			{
				var log2n = log2ns[dn - d - 1];
				var n = 1<<log2n;

				var shift = lowerBits;

				var nSize = 1 << ( totalBits - log2n );
				var upperMask = ~lowerMask << log2n;
				for( var i = 0; i < nSize; ++i )
				{
					var plane = ( i & lowerMask ) | ( i << log2n ) & upperMask;

					FftMultiDPlane( data, shift, plane, log2n );
				}

				lowerBits += log2n;
				lowerMask = ~( ~lowerMask << log2n );
			}
		}

		/// <summary>
		/// Perform an inverse multi-dimensional fast fourier transform.
		/// </summary>
		/// <param name="data">The data to transform.</param>
		/// <param name="log2ns">The log2 order on each dimension.</param>
		public static void InverseMultiFft( this Complex[] data, params byte[] log2ns )
		{
			var nd = log2ns.Length;
			var totalBits = log2ns.Select( b => (int) b ).Sum();

			if( data.Length < ( 1 << totalBits ) )
			{
				throw new ArgumentException( "length/dimension mismatch" );
			}

			var lowerBits = 0;
			var lowerMask = 0;
			for( var d = 0; d < nd; ++d )
			{
				var log2n = log2ns[nd - d - 1];
				var n = 1 << log2n;

				var shift = lowerBits;

				var nSize = 1 << ( totalBits - log2n );
				var upperMask = ~lowerMask << log2n;
				var stride = 1 << lowerBits;
				for( var i = 0; i < nSize; ++i )
				{
					var plane = ( i & lowerMask ) | ( i << log2n ) & upperMask;
					ReverseMultiDPlane( data, plane + stride, n - 1, stride );
					FftMultiDPlane( data, shift, plane, log2n );
				}

				lowerBits += log2n;
				lowerMask = ~( ~lowerMask << log2n );
			}

			var totalLen = data.Length;
			for( var i = 0; i < totalLen; ++i )
			{
				data[i] /= totalLen;
			}
		}

		private static void FftMultiDPlane( Complex[] data, int planeShift, int plane, int log2N )
		{
			var dp = 1 << planeShift;
			var n = 1 << log2N;
			var ipos = plane;
			for( var i = 0; i < n; ++i, ipos += dp )
			{
				var r = ReverseBits( i, log2N );
				if( i >= r ) continue;

				var rpos =  plane + ( r << planeShift );
				(data[ipos], data[rpos]) = (data[rpos], data[ipos]);
			}

			var s = 0;
			var pmax = 1 << ( log2N + planeShift );
			var kmax = plane + pmax;
			var dk = dp + dp;
			for( var dj = dp; dj < pmax; dj = dk, dk += dk )
			{
				var dw = POW2_UNITY_ROOT[s++];
				for( var k0 = plane; k0 < kmax; k0 += dk )
				{
					var k1 = k0 + dj;

					var x0 = data[k0];
					var x1 = data[k1];

					data[k0] = x0 + x1;
					data[k1] = x0 - x1;

					var w = dw;
					var jmax = k0 + dj;
					for( var j0 = k0 + dp; j0 < jmax; w *= dw, j0 += dp )
					{
						var j1 = j0 + dj;

						x0 = data[j0];
						x1 = w * data[j1];

						data[j0] = x0 + x1;
						data[j1] = x0 - x1;
					}
				}
			}
		}

		private static void ReverseMultiDPlane( Complex[] data, int plane, int len, int stride )
		{
			var i = plane;
			var r = plane + ( len - 1 ) * stride;

			while( i < r )
			{
				(data[i], data[r]) = (data[r], data[i]);
				i += stride;
				r -= stride;
			}
		}

		private readonly static Complex[] POW2_UNITY_ROOT = ComputePow2RootsOfUnity();

		private static Complex[] ComputePow2RootsOfUnity()
		{
			const int maxPow2 = 32;
			var result = new Complex[maxPow2];
			result[0] = new Complex( -1, 0 );
			var theta = Math.PI * 0.5;
			for( var i = 1; i < maxPow2; ++i )
			{
				var u = Math.Cos( theta );
				var v = Math.Sin( theta );
				result[i] = new Complex( u, -v );
				theta *= 0.5;
			}
			return result;
		}

		/// <summary>
		/// Reverse the least significant bits in an int value.
		/// </summary>
		/// <param name="n">The int value for bits to reverse.</param>
		/// <param name="bitCount">The number of bits to reverse.</param>
		/// <returns></returns>
		public static int ReverseBits( int n, int bitCount )
		{
			var r = n << ( 32 - bitCount );
			var result = ( BYTE_RBITS[r & 0xFF] << 24 )
				| ( BYTE_RBITS[( r >> 8 ) & 0xFF] << 16 )
				| ( BYTE_RBITS[( r >> 16 ) & 0xFF] << 8 )
				| BYTE_RBITS[( r >> 24 ) & 0xFF];
			return result;
		}

		private static readonly byte[] BYTE_RBITS = CreateByteReversals();
		private static byte[] CreateByteReversals()
		{
			var rbytes = new byte[256];

			for( var i = 0; i < 256; ++i )
			{
				var (fwd, bwd) = (i, 0);
				for( var bit = 0; bit < 8; ++bit )
				{
					bwd <<= 1;
					bwd |= fwd & 1;
					fwd >>= 1;
				}

				rbytes[i] = (byte) bwd;
			}

			return rbytes;
		}

		/// <summary>
		/// Compute floor of log2(n).
		/// </summary>
		/// <param name="value">The input value.</param>
		/// <returns>The bit index of the highest order bit in the value, or -1 if value is zero.</returns>
		public static int FloorLog2( int value )
		{
			var n = value;
			var r = 0;

			if( 0 != ( n & ~0xFFFF ) )
			{
				r += 16;
				n >>= 16;
			}

			if( 0 != ( n & ~0xFF ) )
			{
				r += 8;
				n >>= 8;
			}

			if( 0 != ( n & ~0xF ) )
			{
				r += 4;
				n >>= 4;
			}

			switch( n & 0xF )
			{
				case 0:
					return -1;

				case 1:
					return r;

				case 2:
				case 3:
					return r + 1;

				case 4:
				case 5:
				case 6:
				case 7:
					return r + 2;

				default:
					return r + 3;
			}
		}
	}
}
