using System;
using System.Linq;
using System.Numerics;
using NUnit.Framework;

namespace TrentTobler.Algorithms.FourierTransform.Tests
{
	[TestFixture]
	public class FourierTransformTests
	{
		[Test]
		public void TestDft()
		{
			var dft = new Complex[] { 2, 3, 5, 7, -3, -2, -5, -11 }
				.Dft();

			Assert.AreEqual(
				"-4.00+0.00i | -4.19-26.26i | -1.00-5.00i | 14.19-6.26i | 2.00+0.00i | 14.19+6.26i | -1.00+5.00i | -4.19+26.26i",
				string.Join( " | ", dft.Select( TestExtensions.FormatComplex( "f2" ) ) ),
				"Approximate Values" );
		}

		[Test]
		public void TestDft2DRoundTrip()
		{
			var rand = new Random( 101 );
			var sample = rand.RandomizedComplex2DArray( 8, 16 ).To2D();
			var dft = sample.Dft2D();
			var rt = dft.InverseDft2D();

			var err = (
				from i in Enumerable.Range( 0, sample.GetLength( 0 ) )
				from j in Enumerable.Range( 0, sample.GetLength( 1 ) )
				select ( sample[i, j] - rt[i, j] ).Magnitude
			).Max();

			Assert.IsTrue( err < 1e-10, "round trip error" );
		}

		[Test]
		public void TestDft3DRoundTrip()
		{
			var rand = new Random( 101 );
			var sample = rand.RandomizedComplex3DArray( 8, 8, 8 ).To3D();
			var dft = sample.Dft3D();
			var rt = dft.InverseDft3D();

			var err = (
				from i in Enumerable.Range( 0, sample.GetLength( 0 ) )
				from j in Enumerable.Range( 0, sample.GetLength( 1 ) )
				from k in Enumerable.Range( 0, sample.GetLength( 2 ) )
				select ( sample[i, j, k] - rt[i, j, k] ).Magnitude
			).Max();

			Assert.IsTrue( err < 1e-10, "round trip error" );
		}

		[Test]
		public void TestInverseDft()
		{
			Complex i = Complex.ImaginaryOne;
			var invDft = new Complex[] { -4, -4.19 - 26.26 * i, -1 - 5 * i, 14.19 - 6.26 * i, 2, 14.19 + 6.26 * i, -1 + 5 * i, -4.19 + 26.26 * i }
				.InverseDft();

			var expect = new Complex[] { 2, 3, 5, 7, -3, -2, -5, -11 };
			var err = invDft.GetError( expect );

			Assert.Less( err, 1e-2, string.Concat(
				$"{err} error : ",
				string.Join( " | ", invDft.Select( TestExtensions.FormatComplex( "f2" ) ) ),
				" != ",
				string.Join( " | ", expect.Select( TestExtensions.FormatComplex( "f2" ) ) ) ) );
		}
	}
}
