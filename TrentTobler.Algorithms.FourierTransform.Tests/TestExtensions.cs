using System;
using System.Linq;
using System.Numerics;
using NUnit.Framework;

namespace TrentTobler.Algorithms.FourierTransform.Tests
{
	internal static class TestExtensions
	{
		public static Complex[] RandomizedComplexArray( this Random rand, int len )
		{
			var result = new Complex[len];
			for( var i = 0; i < len; ++i )
			{
				result[i] = new Complex( 4 * rand.NextDouble() - 2, 4 * rand.NextDouble() - 2 );
			}

			return result;
		}

		public static Complex[][] RandomizedComplex2DArray( this Random rand, int len0, int len1 )
		{
			var result = new Complex[len0][];
			for( var i0 = 0; i0 < len0; ++i0 )
			{
				result[i0] = rand.RandomizedComplexArray( len1 );
			}
			return result;
		}

		public static Complex[][][] RandomizedComplex3DArray( this Random rand, int len0, int len1, int len2 )
		{
			var result = new Complex[len0][][];
			for( var i0 = 0; i0 < len0; ++i0 )
			{
				result[i0] = rand.RandomizedComplex2DArray( len1, len2 );
			}
			return result;
		}

		public static Func<Complex, string> FormatComplex( string format )
			=> c =>
			{
				var r = c.Real.ToString( format );
				var i = c.Imaginary.ToString( format );
				return i.StartsWith( "-" ) ? r + i + "i" : r + "+" + i + "i";
			};

		public static Complex[,] To2D( this Complex[][] jaggedArray )
		{
			var len0 = jaggedArray.Length;
			var len1 = jaggedArray[0].Length;

			var result = new Complex[len0, len1];
			for( var i0 = 0; i0 < len0; ++i0 )
			{
				Assert.AreEqual( len1, jaggedArray[i0].Length );
				for( var i1 = 0; i1 < len1; ++i1 )
				{
					result[i0, i1] = jaggedArray[i0][i1];
				}
			}

			return result;
		}

		public static Complex[,,] To3D( this Complex[][][] jaggedArray )
		{
			var len0 = jaggedArray.Length;
			var len1 = jaggedArray[0].Length;
			var len2 = jaggedArray[0][0].Length;

			var result = new Complex[len0, len1, len2];
			for( var i0 = 0; i0 < len0; ++i0 )
			{
				Assert.AreEqual( len1, jaggedArray[i0].Length );
				for( var i1 = 0; i1 < len1; ++i1 )
				{
					Assert.AreEqual( len2, jaggedArray[i0][i1].Length );
					for( var i2 = 0; i2 < len2; ++i2 )
					{
						result[i0, i1, i2] = jaggedArray[i0][i1][i2];
					}
				}
			}

			return result;
		}

		public static Complex[] Flatten2D( this Complex[,] data )
		{
			var result = new Complex[data.Length];
			var len0 = data.GetLength( 0 );
			var len1 = data.GetLength( 1 );
			var p = 0;
			for( var i0 = 0; i0 < len0; ++i0 )
			{
				for( var i1 = 0; i1 < len1; ++i1 )
				{
					result[p] = data[i0, i1];
					++p;
				}
			}
			return result;
		}

		public static Complex[] Flatten3D( this Complex[,,] data )
		{
			var result = new Complex[data.Length];
			var len0 = data.GetLength( 0 );
			var len1 = data.GetLength( 1 );
			var len2 = data.GetLength( 2 );
			var p = 0;
			for( var i0 = 0; i0 < len0; ++i0 )
			{
				for( var i1 = 0; i1 < len1; ++i1 )
				{
					for( var i2 = 0; i2 < len2; ++i2 )
					{
						result[p] = data[i0, i1, i2];
						++p;
					}
				}
			}
			return result;
		}

		public static double GetError( this Complex[] x, Complex[] origX )
			=> origX.Zip( x, ( l, r ) => ( l - r ).Magnitude ).Max();
	}
}
