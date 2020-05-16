using System;
using System.Linq;
using System.Numerics;
using NUnit.Framework;

namespace TrentTobler.Algorithms.FourierTransform.Tests
{
	[TestFixture]
	public class FastFourierTransformTests
	{
		[Test]
		public void TestFft()
		{
			var rand = new Random( 101 );
			for( var i = 0; i < 10; ++i )
			{
				var x = rand.RandomizedComplexArray( 64 );

				var dft = x.Dft();
				x.Fft();

				var err = dft.GetError( x );
				Assert.IsTrue( err < 1e-10, $"error: {err}" );
			}
		}

		[Test]
		public void TestInverseFft()
		{
			var rand = new Random( 101 );
			for( var i = 0; i < 10; ++i )
			{
				var x = rand.RandomizedComplexArray( 64 );

				var invDft = x.InverseDft();
				x.InverseFft();

				var err = invDft.Zip( x, ( l, r ) => ( l - r ).Magnitude ).Max();
				Assert.IsTrue( err < 1e-10, $"error: {err}" );
			}
		}

		[Test]
		public void TestFft2D()
		{
			var rand = new Random( 101 );
			for( var i = 0; i < 10; ++i )
			{
				var data = rand.RandomizedComplex2DArray( 8, 16 );
				var data2D = data.To2D();
				var x = data.SelectMany( row => row ).ToArray();

				Assert.IsTrue(
					x.SequenceEqual( data2D.Flatten2D() ),
					"flatten" );

				var dft = data2D.Dft2D();

				x.MultiFft( 3, 4 );

				var err = dft.Flatten2D().Zip( x, ( l, r ) => ( l - r ).Magnitude ).Max();
				Assert.IsTrue( err < 1e-10, $"error: {err}" );
			}
		}

		[Test]
		public void TestFft3D()
		{
			var rand = new Random( 101 );
			for( var i = 0; i < 5; ++i )
			{
				var data = rand.RandomizedComplex3DArray( 8, 4, 16 );
				var data3D = data.To3D();
				var x = data.SelectMany( row => row ).SelectMany( row => row ).ToArray();

				Assert.IsTrue(
					x.SequenceEqual( data3D.Flatten3D() ),
					"flatten" );

				var dft = data3D.Dft3D();

				x.MultiFft( 3, 2, 4 );

				var err = dft.Flatten3D().Zip( x, ( l, r ) => ( l - r ).Magnitude ).Max();
				Assert.IsTrue( err < 1e-6, $"error: {err}" );
			}
		}

		[Test]
		public void TestInverseMultiFftRoundTrip()
		{
			var rand = new Random( 101 );
			for( var i = 0; i < 5; ++i )
			{
				var data = rand.RandomizedComplex3DArray( 8, 4, 16 );
				var x = data.SelectMany( row => row ).SelectMany( row => row ).ToArray();
				var origX = x.ToArray();

				x.MultiFft( 3, 2, 4 );
				x.InverseMultiFft( 3, 2, 4 );

				var err = x.GetError( origX );
				Assert.IsTrue( err < 1e-6, $"error: {err}" );
			}
		}

		[Test]
		public void TestLargeScaleMultiFft()
		{
			var rand = new Random( 101 );
			var x = rand.RandomizedComplexArray( 64 * 64 * 64 );
			var orig = (Complex[]) x.Clone();

			x.MultiFft( 6, 6, 6 );
			x.InverseMultiFft( 6, 6, 6 );
			var err = x.GetError( orig );
			Assert.IsTrue( err < 1e-6, $"error: {err}" );
		}

		[TestCase( 0x12345678, 32, 0x1E6A2C48 )]
		public void TestReverseBits( int arg, int bitCount, int expect )
		{
			Assert.AreEqual( expect, FastFourierTransform.ReverseBits( arg, bitCount ) );
		}

		[TestCase( 101 )]
		public void TestReverseBitsProperties( int seed )
		{
			var rand = new Random( seed );
			for( var sample = 0; sample < 100; ++sample )
			{
				var bitCount = rand.Next( 1, 31 );
				var bits = rand.Next( 0, 1 << bitCount );
				var rev = FastFourierTransform.ReverseBits( bits, bitCount );
				var revRev = FastFourierTransform.ReverseBits( rev, bitCount );
				Assert.AreEqual( bits, revRev, $"[{bitCount}] {bits:X} => {rev:X} => {revRev:X}" );
			}
		}

		[TestCase( 0, -1 )]
		[TestCase( -1, 31 )]
		[TestCase( 1, 0 )]
		[TestCase( 2, 1 )]
		[TestCase( 4, 2 )]
		[TestCase( 8, 3 )]
		[TestCase( 16, 4 )]
		[TestCase( 32, 5 )]
		[TestCase( 64, 6 )]
		[TestCase( 128, 7 )]
		[TestCase( 256, 8 )]
		[TestCase( 3, 1 )]
		[TestCase( 11, 3 )]
		[TestCase( 250, 7 )]
		public void TestFloorLog2( int n, int expect )
		{
			Assert.AreEqual( expect, FastFourierTransform.FloorLog2( n ) );
		}

		[TestCase( 101 )]
		public void TestFloorLog2Properties( int seed )
		{
			var rand = new Random( seed );
			for( var sample = 0; sample < 100; ++sample )
			{
				var n = rand.Next() ^ ( rand.Next() << 16 );
				if( n == 0 )
				{
					continue;
				}

				var log2n = FastFourierTransform.FloorLog2( n );
				var hibit = 1 << log2n;
				var mask = hibit | ( hibit - 1 );
				Assert.AreEqual( 0, n & ~mask, "higher bits" );
				Assert.AreEqual( hibit, n & hibit, "high bit" );
			}
		}
	}
}
