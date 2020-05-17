using System;
using System.Collections.Generic;
using System.ComponentModel;
using System.Data;
using System.Drawing;
using System.Linq;
using System.Text;
using System.Threading.Tasks;
using System.Windows.Forms;

namespace TrentTobler.Examples.FourierTransform
{
	public partial class MainForm : Form
	{
		public MainForm()
		{
			InitializeComponent();
			this.SetStyle(
				ControlStyles.AllPaintingInWmPaint
				| ControlStyles.Opaque
				| ControlStyles.OptimizedDoubleBuffer
				| ControlStyles.ResizeRedraw
				| ControlStyles.UserPaint,
				 true );
		}
	}
}
