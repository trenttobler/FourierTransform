using System;
using System.Numerics;

namespace TrentTobler.Algorithms.FourierTransform
{
	public static class DiscreteFourierTransform
	{
		public static Complex[] Dft( this Complex[] x )
		{
			var len = x.Length;

			var y = new Complex[len];

			for( var k = 0; k < len; ++k )
			{
				var dw = 2 * Math.PI * k / len;

				var sum = x[0];
				for( var n = 1; n < len; ++n )
				{
					var w = dw * n;
					sum += x[n] * new Complex( Math.Cos( w ), -Math.Sin( w ) );
				}

				y[k] = sum;
			}

			return y;
		}

		public static Complex[] InverseDft( this Complex[] x )
		{
			var len = x.Length;

			var y = new Complex[len];

			for( var k = 0; k < len; ++k )
			{
				var dw = 2 * Math.PI * k / len;

				var sum = x[0];
				for( var n = 1; n < len; ++n )
				{
					var w = dw * n;
					sum += x[n] * new Complex( Math.Cos( w ), Math.Sin( w ) );
				}

				y[k] = sum / len;
			}

			return y;
		}

		public static Complex[,] Dft2D( this Complex[,] x )
		{
			var len0 = x.GetLength( 0 );
			var len1 = x.GetLength( 1 );

			var y = new Complex[len0, len1];

			for( var k0 = 0; k0 < len0; ++k0 )
			{
				var dw0 = 2 * Math.PI * k0 / len0;

				for( var k1 = 0; k1 < len1; ++k1 )
				{
					var dw1 = 2 * Math.PI * k1 / len1;

					var sum = Complex.Zero;
					for( var n0 = 0; n0 < len0; ++n0 )
					{
						var w0 = dw0 * n0;

						for( var n1 = 0; n1 < len1; ++n1 )
						{
							var w1 = dw1 * n1;
							var wct = new Complex( Math.Cos( w0 + w1 ), -Math.Sin( w0 + w1 ) );
							sum += x[n0, n1] * wct;
						}
					}

					y[k0, k1] = sum;
				}
			}

			return y;
		}

		public static Complex[,] InverseDft2D( this Complex[,] x )
		{
			var len0 = x.GetLength( 0 );
			var len1 = x.GetLength( 1 );
			var lenT = len0 * len1;

			var y = new Complex[len0, len1];

			for( var k0 = 0; k0 < len0; ++k0 )
			{
				var dw0 = 2 * Math.PI * k0 / len0;

				for( var k1 = 0; k1 < len1; ++k1 )
				{
					var sum = Complex.Zero;
					var dw1 = 2 * Math.PI * k1 / len1;

					for( var n0 = 0; n0 < len0; ++n0 )
					{
						var w0 = dw0 * n0;

						for( var n1 = 0; n1 < len1; ++n1 )
						{
							var w1 = dw1 * n1;
							var wct = new Complex( Math.Cos( w0 + w1 ), Math.Sin( w0 + w1 ) );
							sum += x[n0, n1] * wct;
						}
					}

					y[k0, k1] = sum / lenT;
				}
			}

			return y;
		}

		public static Complex[,,] Dft3D( this Complex[,,] x )
		{
			var len0 = x.GetLength( 0 );
			var len1 = x.GetLength( 1 );
			var len2 = x.GetLength( 2 );

			var y = new Complex[len0, len1, len2];

			for( var k0 = 0; k0 < len0; ++k0 )
			{
				var dw0 = 2 * Math.PI * k0 / len0;

				for( var k1 = 0; k1 < len1; ++k1 )
				{
					var dw1 = 2 * Math.PI * k1 / len1;

					for( var k2 = 0; k2 < len2; ++k2 )
					{
						var dw2 = 2 * Math.PI * k2 / len2;

						var sum = Complex.Zero;

						for( var n0 = 0; n0 < len0; ++n0 )
						{
							var w0 = dw0 * n0;

							for( var n1 = 0; n1 < len1; ++n1 )
							{
								var w1 = dw1 * n1;

								for( var n2 = 0; n2 < len2; ++n2 )
								{
									var w2 = dw2 * n2;

									var wct = new Complex( Math.Cos( w0 + w1 + w2 ), -Math.Sin( w0 + w1 + w2 ) );
									sum += x[n0, n1, n2] * wct;
								}
							}
						}

						y[k0, k1, k2] = sum;
					}
				}
			}

			return y;
		}

		public static Complex[,,] InverseDft3D( this Complex[,,] x )
		{
			var len0 = x.GetLength( 0 );
			var len1 = x.GetLength( 1 );
			var len2 = x.GetLength( 2 );

			var lenT = len0 * len1 * len2;

			var y = new Complex[len0, len1, len2];

			for( var k0 = 0; k0 < len0; ++k0 )
			{
				var dw0 = 2 * Math.PI * k0 / len0;

				for( var k1 = 0; k1 < len1; ++k1 )
				{
					var dw1 = 2 * Math.PI * k1 / len1;

					for( var k2 = 0; k2 < len2; ++k2 )
					{
						var dw2 = 2 * Math.PI * k2 / len2;

						var sum = Complex.Zero;

						for( var n0 = 0; n0 < len0; ++n0 )
						{
							var w0 = dw0 * n0;

							for( var n1 = 0; n1 < len1; ++n1 )
							{
								var w1 = dw1 * n1;

								for( var n2 = 0; n2 < len2; ++n2 )
								{
									var w2 = dw2 * n2;

									var wct = new Complex( Math.Cos( w0 + w1 + w2 ), Math.Sin( w0 + w1 + w2 ) );
									sum += x[n0, n1, n2] * wct;
								}
							}
						}

						y[k0, k1, k2] = sum / lenT;
					}
				}
			}

			return y;
		}
	}
}
