using System;
using System.Linq;
using System.Windows.Forms;
using System.Drawing;
using System.Drawing.Imaging;
using System.Numerics;
using System.Runtime.InteropServices;
using TrentTobler.Algorithms.FourierTransform;
using System.Diagnostics;

namespace TrentTobler.Examples.FourierTransform
{
	static class Program
	{
		[STAThread]
		static void Main()
		{
			Application.SetHighDpiMode( HighDpiMode.SystemAware );
			Application.EnableVisualStyles();
			Application.SetCompatibleTextRenderingDefault( false );

			var dx = 256;
			var dy = 256;
			var dz = 256;

			var data = GenerateCloudFft( dx, dy, dz );

			var mainForm = new MainForm();
			var timer = new Timer
			{
				Interval = 50
			};

			var bitmap = new Bitmap( dx, dy, PixelFormat.Format24bppRgb );
			byte[] bitmapBytes = null;
			var framez = 0;

			GenerateCloudBitmap( framez, dx, dy, data, bitmap, ref bitmapBytes );

			timer.Tick += ( sender, e ) =>
			{
				if( ++framez >= dz )
				{
					framez = 0;
				}

				GenerateCloudBitmap( framez, dx, dy, data, bitmap, ref bitmapBytes );

				mainForm.Invalidate();
			};

			mainForm.Load += ( sender, e ) => timer.Start();
			mainForm.Paint += ( sender, e ) => RenderCloudBitmap( e, mainForm, bitmap );

			Application.Run( mainForm );
		}

		private static float[] GenerateCloudFft( int dx, int dy, int dz )
		{
			var log2X = (byte) FastFourierTransform.FloorLog2( dx );
			var log2Y = (byte) FastFourierTransform.FloorLog2( dy );
			var log2Z = (byte) FastFourierTransform.FloorLog2( dz );

			Trace.WriteLine( $"Generate random {dx} x {dy} x {dz} data..." );
			var data = new Complex[dx * dy * dz];
			var rand = new Random( 101 );
			for( var i = 0; i < data.Length; ++i )
			{
				data[i] = rand.NextDouble();
			}

			Trace.WriteLine( $"Compute multi-dimensional FFT..." );
			data.MultiFft( log2X, log2Y, log2Z );

			Trace.WriteLine( $"Apply 1/f^2 scaling..." );
			var p = 0;
			for( var x = 0; x < dx; ++x )
			{
				var fx = Math.Min( x, dx - x );
				for( var y = 0; y < dy; ++y )
				{
					var fy = Math.Min( y, dy - y );

					for( var z = 0; z < dz; ++z )
					{
						var fz = Math.Min( z, dz - z );

						var ff = fx * fx + fy * fy + fz * fz;
						if( ff == 0 )
						{
							data[p++] = 0;
						}
						else
						{
							data[p++] /= ff;
						}
					}
				}
			}

			Trace.WriteLine( "Compute inverse multi-dimensional FFT..." );
			data.InverseMultiFft( log2X, log2Y, log2Z  );

			Trace.WriteLine( "Convert complex real component to normalized float..." );
			var minc = data.Min( x => x.Real );
			var maxc = data.Max( x => x.Real );
			var result = new float[data.Length];
			for( var i = 0; i < data.Length; ++i )
			{
				result[i] = (float)( ( data[i].Real - minc ) / ( maxc - minc ) );
			}

			Trace.WriteLine( $"max( abs( R ) ) = {data.Max( c => Math.Abs( c.Real ) )}" );
			Trace.WriteLine( $"max( abs( I ) ) = {data.Max( c => Math.Abs( c.Imaginary ) )}" );

			return result;
		}

		private static void GenerateCloudBitmap(
			int z,
			int dx, int dy,
			float[] data,
			Bitmap bitmap,
			ref byte[] bitmapBytes )
		{
			var lockedBits = bitmap.LockBits( new Rectangle( 0, 0, dx, dy ), ImageLockMode.WriteOnly, PixelFormat.Format24bppRgb );

			var pixelDataLen = Math.Abs( lockedBits.Stride ) * dy;
			if( bitmapBytes == null || bitmapBytes.Length < pixelDataLen )
			{
				bitmapBytes = new byte[pixelDataLen];
			}

			var dplane = dx * dy;
			var inp = z * dplane;
			var outp = 0;
			for( var i = 0; i < dplane; ++i )
			{
				var c = data[inp++];

				bitmapBytes[outp++] = GetPixelByte( 128 + c * 127 ); // blue
				bitmapBytes[outp++] = GetPixelByte( c * c * 256 ); // green
				bitmapBytes[outp++] = GetPixelByte( c * c * 256 ); // red
			}

			var ptr0 = lockedBits.Stride >= 0 ? lockedBits.Scan0 : lockedBits.Scan0 + lockedBits.Stride * ( dy - 1 );
			Marshal.Copy( bitmapBytes, 0, ptr0, pixelDataLen );

			bitmap.UnlockBits( lockedBits );
		}
		private static void RenderCloudBitmap( PaintEventArgs e, MainForm mainForm, Bitmap bitmap )
		{
			e.Graphics.Clear( Color.Black );

			var d = Math.Min( mainForm.ClientSize.Width, mainForm.ClientSize.Height );
			var x0 = ( mainForm.ClientSize.Width - d ) / 2;
			var y0 = ( mainForm.ClientSize.Height - d ) / 2;
			e.Graphics.DrawImage( bitmap, new Rectangle( x0, y0, d, d ) );
		}

		private static byte GetPixelByte( double c )
			=> (byte) Math.Max( 0, Math.Min( 255, (int) c ) );
	}
}
