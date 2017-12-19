package plot;
import java.awt.Color;

import javax.swing.JFrame;

import java.awt.BasicStroke; 

import org.jfree.chart.ChartPanel; 
import org.jfree.chart.JFreeChart; 
import org.jfree.data.xy.XYDataset; 
import org.jfree.data.xy.XYSeries; 
import org.jfree.ui.ApplicationFrame; 
import org.jfree.ui.RefineryUtilities; 
import org.jfree.chart.plot.XYPlot; 
import org.jfree.chart.ChartFactory; 
import org.jfree.chart.plot.PlotOrientation; 
import org.jfree.data.xy.XYSeriesCollection; 
import org.jfree.chart.renderer.xy.XYLineAndShapeRenderer;
public class XyChart extends JFrame {
	public static XYDataset dataSet; 
	public XyChart(String applicationTitle, String chartTitle ) {
		super(applicationTitle);
		// TODO Auto-generated constructor stub
		JFreeChart xylineChart = ChartFactory.createXYLineChart(
		         chartTitle ,
		         "Feature Index" ,
		         "Feature" ,
		         dataSet ,
		         PlotOrientation.VERTICAL ,
		         true , true , false);
		         
		      ChartPanel chartPanel = new ChartPanel( xylineChart );
		      chartPanel.setPreferredSize( new java.awt.Dimension( 560 , 367 ) );
		      final XYPlot plot = xylineChart.getXYPlot( );
		      XYLineAndShapeRenderer renderer = new XYLineAndShapeRenderer( );
		      renderer.setSeriesPaint( 0 , Color.RED );
		      renderer.setSeriesPaint( 1 , Color.BLUE );
		      //renderer.setSeriesPaint( 2 , Color.YELLOW );
		      renderer.setSeriesStroke( 0 , new BasicStroke( 1.0f ) );
		      renderer.setSeriesStroke( 1 , new BasicStroke( 1.0f ) );
		      //renderer.setSeriesStroke( 2 , new BasicStroke( 1.0f ) );
		      plot.setRenderer( renderer ); 
		      setContentPane( chartPanel ); 
		      this.setDefaultCloseOperation(JFrame.DISPOSE_ON_CLOSE);
	}
	

}
