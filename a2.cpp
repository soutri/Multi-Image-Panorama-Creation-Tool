// B657 assignment 2 skeleton code
//
// Compile with: "make"
//
// See assignment handout for command line and project specifications.


//Link to the header file
#include "CImg.h"
#include <iostream>
#include <stdlib.h>
#include <string>
#include <vector>
#include <Sift.h>


//Use the cimg namespace to access functions easily
using namespace cimg_library;
using namespace std;
struct Point {
  double x,y;
};

struct PointPair {
  double x1,x2,y1,y2;
};
struct MatchPointInfo {
  double x,y,dist;
};
struct homogdist {
  CImg<double> homog;
  int match_count;
  vector<PointPair> good_matches;
};
struct trans_images_struct
{
	CImg<double> image;
};

void draw_descriptor_image(CImg<double> image, const vector<SiftDescriptor> descriptors, const char *filename)
{
  for(unsigned int i=0; i < descriptors.size(); i++)
    {
      int tx1 = 0, ty1 = 0, tx2 = 0, ty2 = 0;
      double color_point[] = {255.0, 255.0, 0};
      for(int x=-2; x<3; x++)
	for(int y=-2; y<3; y++)
	  if(x==0 || y==0)
	    for(int c=0; c<3; c++){
	      //Find if coordinates are in workspace to draw crosshair
	      tx1 = (descriptors[i].col + y - 1);
	      ty1 = (descriptors[i].row + x - 1);
	      if (tx1 >= 0 && tx1 < image.width() && ty1 >= 0 && ty1 < image.height())
		image( tx1, ty1, 0, c) = color_point[c];				
	    }
    }
  image.get_normalize(0,255).save(filename);
}

CImg<double> apply_transformation(CImg<double> inputimage,CImg<double> transformation)
{
	CImg<double> output_image(inputimage.width(),inputimage.height(),inputimage.depth(),inputimage.spectrum(),0);
	CImg<double> transformation_inv=transformation.invert(true);
	const int width = inputimage.width();
	const int height = inputimage.height();
	
	for(int x=0;x<width;x++)
	{
		for(int y=0;y<height;y++)
		{
			double new_x =((transformation_inv(0,0)*x) + (transformation_inv(0,1)*y) + transformation_inv(0,2));
			double new_y = ((transformation_inv(1,0)*x)+ (transformation_inv(1,1)*y) + transformation_inv(1,2)) ;
            double w =((transformation_inv(2,0)*x)+ (transformation_inv(2,1)*y)+transformation_inv(2,2));
			//cout<<new_x<<"\t"<<new_y<<"\t"<<w<<"\n";
			if(w!=0)
			{
				new_x=(new_x/w);
				new_y=(new_y/w);
			}
			//cout<<new_x<<"\t"<<new_y<<"\n";
			if(new_x>=0 && new_x<width && new_y>=0 && new_y<height)
			{
				output_image.atXY(x,y,0)=inputimage.atXY(new_x,new_y,0);
				output_image.atXY(x,y,1)=inputimage.atXY(new_x,new_y,1);
				output_image.atXY(x,y,2)=inputimage.atXY(new_x,new_y,2);
			}
		}
		
	}

	return output_image;
}

CImg<double> apply_transformation2(CImg<double> output_image,CImg<double> inputimage,CImg<double> transformation,double offset1_x,double offset1_y,double offset2,bool do_check)
{
	CImg<double> transformation_inv=transformation.invert(true);
	const int width = inputimage.width();
	const int height = inputimage.height();
	
	for(int x=0;x<width;x++)
	{
		for(int y=0;y<height;y++)
		{
			double new_x =((transformation_inv(0,0)*x) + (transformation_inv(0,1)*y) + transformation_inv(0,2));
			double new_y = ((transformation_inv(1,0)*x)+ (transformation_inv(1,1)*y) + transformation_inv(1,2)) ;
            double w =((transformation_inv(2,0)*x)+ (transformation_inv(2,1)*y)+transformation_inv(2,2));
			//cout<<new_x<<"\t"<<new_y<<"\t"<<w<<"\n";
			if(w!=0)
			{
				new_x=(new_x/w);
				new_y=(new_y/w);
			}
			//cout<<new_x<<"\t"<<new_y<<"\n";
			if(new_x+offset1_x>=0 && new_x+offset1_x<width && new_y+offset1_y>=0 && new_y+offset1_y<height)
			{
				int r=output_image.atXY(x+offset2,y,0);
				int g=output_image.atXY(x+offset2,y,1);
				int b=output_image.atXY(x+offset2,y,2);
				bool check=(r==0&&b==0&&g==0);
				if(true&&do_check)
				{
					output_image.atXY(x+offset2,y,0)=inputimage.atXY(new_x+offset1_x,new_y+offset1_y,0);
					output_image.atXY(x+offset2,y,1)=inputimage.atXY(new_x+offset1_x,new_y+offset1_y,1);
					output_image.atXY(x+offset2,y,2)=inputimage.atXY(new_x+offset1_x,new_y+offset1_y,2);
				}
			}
		}
		
	}

	return output_image;
}

struct myclass {
  bool operator() (MatchPointInfo i,MatchPointInfo j) { return (i.dist<j.dist);}
} myobject;

vector<PointPair> find_matches_in_images(CImg<double> &image1,CImg<double> &image2,double thresh)
{
	CImg<double> image1_gray = image1.get_RGBtoHSI().get_channel(2);   
	vector<SiftDescriptor> image1_descriptors = Sift::compute_sift(image1_gray);
	CImg<double> image2_gray = image2.get_RGBtoHSI().get_channel(2);   
	vector<SiftDescriptor> image2_descriptors = Sift::compute_sift(image2_gray);
	
	int size1=image1_descriptors.size();
	int size2=image2_descriptors.size();
	vector<PointPair> sift_matches;
	for(int i=0;i<size1;i++)
	{
		vector<Point> matches;
		vector<MatchPointInfo> distances;
		for(int j=0;j<size2;j++)
		{
			double dist=0;
			//looping through the individual descriptor
			for(int k=0;k<128;k++)
				dist+= pow(image1_descriptors[i].descriptor[k] - image2_descriptors[j].descriptor[k],2);
			dist=sqrt(dist);
			MatchPointInfo matching_point={image2_descriptors[j].col,image2_descriptors[j].row,dist};
			distances.push_back(matching_point);
		}
		sort(distances.begin(),distances.end(),myobject);
		MatchPointInfo best_match_1=distances[0];
		MatchPointInfo best_match_2=distances[1];
		if(best_match_1.dist/best_match_2.dist<thresh)
		{
			struct PointPair temp = {image1_descriptors[i].col,best_match_1.x,image1_descriptors[i].row,best_match_1.y};
			sift_matches.push_back(temp);
		}
	}
	cout<<"Total sift Matches:"<<sift_matches.size()<<"\n";
	//for(int l=0;l<sift_matches.size();l++)
	//	cout<<sift_matches[l].x1<<"\t"<<sift_matches[l].y1<<"\t"<<sift_matches[l].x2<<"\t"<<sift_matches[l].y2<<"\n";
	return sift_matches;
}
CImg<double> draw_lines_sift_matches(CImg<double> desc_img1,CImg<double> desc_img2,vector<PointPair> sift_matches)
{
	  CImg<double> output_image(desc_img1.width()+desc_img2.width(),desc_img1.height(),desc_img1.depth(),desc_img1.spectrum(),0);
	  for(int x=0; x<=desc_img1.width();x++)
	  {
		  for(int y=0;y<=desc_img1.height();y++)
		  {
			  output_image.atXY(x,y,0)=desc_img1.atXY(x,y,0);
			  output_image.atXY(x,y,1)=desc_img1.atXY(x,y,1);
			  output_image.atXY(x,y,2)=desc_img1.atXY(x,y,2);
		  }
	  }
	  for(int x=0; x<=desc_img2.width();x++)
	  {
		  for(int y=0;y<=desc_img2.height();y++)
		  {
			  output_image.atXY(x+desc_img1.width(),y,0)=desc_img2.atXY(x,y,0);
			  output_image.atXY(x+desc_img1.width(),y,1)=desc_img2.atXY(x,y,1);
			  output_image.atXY(x+desc_img1.width(),y,2)=desc_img2.atXY(x,y,2);
		  }
	  }
	  output_image.save("joined.png");
	 
	  CImg<double> new_img("joined.png");
	  const unsigned char color[] = {250,250,210};
	  for(int l=0;l<sift_matches.size();l++)
		new_img.draw_line(sift_matches[l].x1,sift_matches[l].y1,sift_matches[l].x2+desc_img1.width(),sift_matches[l].y2,color); 
		 
	  return new_img;
}
//Referred http://www.csc.kth.se/~perrose/files/pose-init-model/node17.html for the linear equation solving for
//finding the homography_matrix
CImg<double> solve_homography(vector<Point> image1,vector<Point> image2)
{
	CImg<double> eq_mat(8,8,1,1);
	
	double x_1=image1[0].x,y_1=image1[0].y,x1=image2[0].x,y1=image2[0].y;
	double x_2=image1[1].x,y_2=image1[1].y,x2=image2[1].x,y2=image2[1].y;
	double x_3=image1[2].x,y_3=image1[2].y,x3=image2[2].x,y3=image2[2].y;
	double x_4=image1[3].x,y_4=image1[3].y,x4=image2[3].x,y4=image2[3].y;
	
    eq_mat(0,0) = x1;	eq_mat(1,0) = y1;	eq_mat(2,0) = 1;	eq_mat(3,0) = 0; 	eq_mat(4,0) = 0; 	eq_mat(5,0) = 0;	eq_mat(6,0) = (-x1*x_1);	eq_mat(7,0) = (x_1*(-y1));
    eq_mat(0,1) = 0; 	eq_mat(1,1) = 0; 	eq_mat(2,1) = 0;	eq_mat(3,1) = x1;	eq_mat(4,1) = y1;	eq_mat(5,1) = 1;	eq_mat(6,1) = (-x1*y_1);	eq_mat(7,1) = (y_1*(-y1));
    eq_mat(0,2) = x2;	eq_mat(1,2) = y2;	eq_mat(2,2) = 1;	eq_mat(3,2) = 0; 	eq_mat(4,2) = 0; 	eq_mat(5,2) = 0;	eq_mat(6,2) = (-x2*x_2);	eq_mat(7,2) = (x_2*(-y2));
    eq_mat(0,3) = 0; 	eq_mat(1,3) = 0; 	eq_mat(2,3) = 0;	eq_mat(3,3) = x2;	eq_mat(4,3) = y2;	eq_mat(5,3) = 1;	eq_mat(6,3) = (-x2*y_2);	eq_mat(7,3) = (y_2*(-y2));
    eq_mat(0,4) = x3;	eq_mat(1,4) = y3;	eq_mat(2,4) = 1;	eq_mat(3,4) = 0; 	eq_mat(4,4) = 0; 	eq_mat(5,4) = 0;	eq_mat(6,4) = (-x3*x_3);	eq_mat(7,4) = (x_3*(-y3));
    eq_mat(0,5) = 0; 	eq_mat(1,5) = 0; 	eq_mat(2,5) = 0;	eq_mat(3,5) = x3;	eq_mat(4,5) = y3;	eq_mat(5,5) = 1;	eq_mat(6,5) = (-x3*y_3);	eq_mat(7,5) = (y_3*(-y3));
    eq_mat(0,6) = x4;	eq_mat(1,6) = y4;	eq_mat(2,6) = 1;	eq_mat(3,6) = 0; 	eq_mat(4,6) = 0; 	eq_mat(5,6) = 0;	eq_mat(6,6) = (-x4*x_4);	eq_mat(7,6) = (x_4*(-y4));
    eq_mat(0,7) = 0; 	eq_mat(1,7) = 0; 	eq_mat(2,7) = 0;	eq_mat(3,7) = x4;	eq_mat(4,7) = y4;	eq_mat(5,7) = 1;	eq_mat(6,7) = (-x4*y_4);	eq_mat(7,7) = (y_4*(-y4));

    CImg<double> inverse_mat = eq_mat.invert(true);
	CImg<double> points(1,8,1,1,x_1,y_1,x_2,y_2,x_3,y_3,x_4,y_4);

    inverse_mat*=points;

    CImg<double> homography_matrix(3,3,1,1,inverse_mat(0,0),inverse_mat(1,0),inverse_mat(2,0),inverse_mat(3,0),inverse_mat(4,0),inverse_mat(5,0),inverse_mat(6,0),inverse_mat(7,0),1.0);

	return homography_matrix.get_transpose();
}


CImg<double> create_billboard(CImg<double> inputimage,vector<Point> billboard_corners,CImg<double> poster)
{
	CImg<double> output_image(inputimage.width(),inputimage.height(),inputimage.depth(),inputimage.spectrum(),0);

	vector<Point> poster_corners;
	poster_corners.push_back({0,0});
	poster_corners.push_back({poster.width(),0});
	poster_corners.push_back({0,poster.height()});
	poster_corners.push_back({poster.width(),poster.height()});
	CImg<double> homography_matrix= solve_homography(billboard_corners,poster_corners);
	
	int min_x=billboard_corners[0].x;
	int max_x=billboard_corners[3].x;
	int min_y=billboard_corners[0].y<billboard_corners[1].y?billboard_corners[0].y:billboard_corners[1].y;
	int max_y=billboard_corners[2].y>billboard_corners[3].y?billboard_corners[2].y:billboard_corners[3].y;
	//output_image=inputimage;
	CImg<double> transformation_inv=homography_matrix.invert(true);
	for(int x=0;x<inputimage.width();x++)
	{
		for(int y=0;y<inputimage.height();y++)
		{
			double new_x =((transformation_inv(0,0)*x) + (transformation_inv(0,1)*y) + transformation_inv(0,2));
			double new_y = ((transformation_inv(1,0)*x)+ (transformation_inv(1,1)*y) + transformation_inv(1,2)) ;
            double w =((transformation_inv(2,0)*x)+ (transformation_inv(2,1)*y)+transformation_inv(2,2));
			//cout<<new_x<<"\t"<<new_y<<"\t"<<w<<"\n";
			int r=inputimage.atXY(x,y,0);
			int g=inputimage.atXY(x,y,1);
			int b=inputimage.atXY(x,y,2);
			bool check=(r+g+b>=750);
			if(check && x>=(min_x-30) && x<(max_x+30) && y>=(min_y) && y<(max_y))
			{
			if(w!=0)
			{
				new_x=(new_x/w);
				new_y=(new_y/w);
			}
			//cout<<new_x<<"\t"<<new_y<<"\n";

				output_image.atXY(x,y,0)=poster.atXY(new_x,new_y,0);
				output_image.atXY(x,y,1)=poster.atXY(new_x,new_y,1);
				output_image.atXY(x,y,2)=poster.atXY(new_x,new_y,2);
			}
			else
			{
				output_image.atXY(x,y,0)=inputimage.atXY(x,y,0);
				output_image.atXY(x,y,1)=inputimage.atXY(x,y,1);
				output_image.atXY(x,y,2)=inputimage.atXY(x,y,2);
		
			}
		}
		
	}
	return output_image;

}
struct myclass1 {
  bool operator() (homogdist i,homogdist j) { return (i.match_count<j.match_count);}
} myobject1;
struct homogdist ransac(vector<PointPair> matches,double dist)
{
	int total_matches=matches.size();
	int index;
	int iterations=500;
	vector<PointPair> good_matches;
	vector<homogdist> homogs;
	//CImg<double> homography_matrix(3,3,1,1);
	//for(int i=0;i<200;i++)
	//	homogs.push_back({homography_matrix,0});
	for(int iter=0;iter<iterations;iter++)
	{
		//cout<<"Iteration:"<<iter<<"\n";
		vector<Point> img1_points;
		vector<Point> img2_points;
		for(int i=0;i<4;i++)
		{
			index=rand()%total_matches;
			struct PointPair temp = matches[index];
			img1_points.push_back({temp.x1,temp.y1});
			img2_points.push_back({temp.x2,temp.y2});
			//cout<<temp.x1<<"\t"<<temp.y1<<"\t"<<temp.x2<<"\t"<<temp.y2<<"\n";

		}
		CImg<double> homography_matrix= solve_homography(img1_points,img2_points);
		CImg<double> transformation_inv=homography_matrix.invert(true);
		int inliers=0;
		homogs.push_back({homography_matrix,0});
		//cout<<"Matches:"<<"\n";
		for(int l=0;l<matches.size();l++)
		{
			struct PointPair temp = matches[l];
			double x=temp.x1,y=temp.y1,target_x=temp.x2,target_y=temp.y2;
			double new_x =((transformation_inv(0,0)*x) + (transformation_inv(0,1)*y) + transformation_inv(0,2));
			double new_y = ((transformation_inv(1,0)*x)+ (transformation_inv(1,1)*y) + transformation_inv(1,2)) ;
			double w =((transformation_inv(2,0)*x)+ (transformation_inv(2,1)*y)+transformation_inv(2,2));
			if(w!=0)
			{
				new_x=(new_x/w);
				new_y=(new_y/w);
			}
			//cout<<new_x<<"\t"<<target_x<<"\t"<<new_y<<"\t"<<target_y<<"\n";
			dist= pow(target_x - new_x,2)+pow(target_y - new_y,2);
			if(dist<=4)
			{
				//cout<<new_x<<"\t"<<target_x<<"\t"<<new_y<<"\t"<<target_y<<"\n";
				homogs[iter].match_count+=1;
				homogs[iter].good_matches.push_back(matches[l]);
			}
				
		}
	}
	//for(int i=0;i<homogs.size();i++)
	//	cout<<homogs[i].match_count<<"\t";
	
	sort(homogs.begin(),homogs.end(),myobject1);
	
	CImg<double> best_homog=homogs[iterations-1].homog;
	
	int best_count=homogs[iterations-1].match_count;
	
	//cout<<"Best match count:"<<best_count<<"\n";
	
	good_matches=homogs[iterations-1].good_matches;
	
	return homogs[iterations-1];
}


vector<trans_images_struct> transform_images(CImg<double> image1,CImg<double> image2,CImg<double> image3)
{
	double thresh=0.8;
	vector<trans_images_struct> trans_images;
	vector<PointPair> matches_1=find_matches_in_images(image3,image2,thresh);
	homogdist ransac_output_1=ransac(matches_1,5.0);
	//CImg<double> good_match_image=draw_lines_sift_matches(image2,image3,ransac_output.good_matches);
	//good_match_image.save("stitch_image_match.png");
	CImg<double> homog_1=ransac_output_1.homog;
	CImg<double> final_image(image1.width()*3,image1.height()*2,image1.depth(),image1.spectrum(),0);
	CImg<double> translation(3,3,1,1, 1,0,final_image.width()/2,0,1,0,0,0,1);	
	//homog_1*=translation;
	CImg<double> transformation_inv=homog_1.invert(true);
	int x=0,y=0;
	double new_x =((transformation_inv(0,0)*x) + (transformation_inv(0,1)*y) + transformation_inv(0,2));
	double new_y = ((transformation_inv(1,0)*x)+ (transformation_inv(1,1)*y) + transformation_inv(1,2)) ;
    double w =((transformation_inv(2,0)*x)+ (transformation_inv(2,1)*y)+transformation_inv(2,2));
	//cout<<new_x<<"\t"<<new_y<<"\t"<<w<<"\n";
	if(w!=0)
	{
		new_x=(new_x/w);
		new_y=(new_y/w);
	}
	//cout<<new_x<<"\t"<<new_y<<"\n";
	double offset1_x=-1*new_x;
	double offset1_y=0;	
	CImg<double> warped_image1=apply_transformation2(final_image,image3,ransac_output_1.homog,offset1_x,0,0.95*image1.width(),true);
				
	trans_images.push_back({warped_image1});
	vector<PointPair> matches_3=find_matches_in_images(image2,warped_image1,thresh);
	homogdist ransac_output_3=ransac(matches_3,5.0);
	CImg<double> homog_3=ransac_output_3.homog;
	//homog_3*=translation;
	transformation_inv=homog_3.invert(true);
	
	x=0,y=0;
	new_x =((transformation_inv(0,0)*x) + (transformation_inv(0,1)*y) + transformation_inv(0,2));
	new_y = ((transformation_inv(1,0)*x)+ (transformation_inv(1,1)*y) + transformation_inv(1,2)) ;
    w =((transformation_inv(2,0)*x)+ (transformation_inv(2,1)*y)+transformation_inv(2,2));
	//cout<<new_x<<"\t"<<new_y<<"\t"<<w<<"\n";
	if(w!=0)
	{
		new_x=(new_x/w);
		new_y=(new_y/w);
	}
	//cout<<new_x<<"\t"<<new_y<<"\n";
	double offset2_x=-1*new_x;
	double offset2_y=0;	
	CImg<double> warped_image3=apply_transformation2(warped_image1,image2,ransac_output_3.homog,offset2_x,0,0.63*image1.width(),true);
	trans_images.push_back({warped_image3});
	
	
	vector<PointPair> matches_2=find_matches_in_images(image1,warped_image3,thresh);
	homogdist ransac_output_2=ransac(matches_2,5.0);
	CImg<double> homog_2=ransac_output_2.homog;
	transformation_inv=homog_2.invert(true);
	
	x=0,y=0;
	new_x =((transformation_inv(0,0)*x) + (transformation_inv(0,1)*y) + transformation_inv(0,2));
	new_y = ((transformation_inv(1,0)*x)+ (transformation_inv(1,1)*y) + transformation_inv(1,2)) ;
    w =((transformation_inv(2,0)*x)+ (transformation_inv(2,1)*y)+transformation_inv(2,2));
	//cout<<new_x<<"\t"<<new_y<<"\t"<<w<<"\n";
	if(w!=0)
	{
		new_x=(new_x/w);
		new_y=(new_y/w);
	}
	//cout<<new_x<<"\t"<<new_y<<"\n";
	CImg<double> warped_image2=apply_transformation2(warped_image3,image1,ransac_output_2.homog,0,0,0,true);
	trans_images.push_back({warped_image2});

	return trans_images;

}
CImg<double> Gaussian_Pyramid(CImg<double> prev_layer, CImg<double> kernel)
{
	CImg<double> prev_layer_conv = prev_layer.get_convolve(kernel);
	CImg<double> current_layer = prev_layer_conv.resize_halfXY();	
	return current_layer;
}

CImg<double> Difference(CImg<double> Gaussian_layer_i, CImg<double> Laplacian_Gaussian_iplus1)
{
	int rows = Gaussian_layer_i.height();
	int cols = Gaussian_layer_i.width();
	//Gaussian_layer_i = Gaussian_layer_i.normalize(0,255);
	//Laplacian_Gaussian_iplus1 = Laplacian_Gaussian_iplus1.normalize(0,255);
	CImg<double> subtracted_img(rows, cols,Gaussian_layer_i.depth(), Gaussian_layer_i.spectrum());
	for(int c = 0; c<3; c++)
		for(int i = 0; i<rows; i++)
			for(int j = 0; j<cols; j++)
				subtracted_img.atXY(i,j,c) = Gaussian_layer_i.atXY(i,j,c) - Laplacian_Gaussian_iplus1.atXY(i,j,c);
	return subtracted_img;
}

CImg<double> Blended_Laplacian(CImg<double> mask_pyramid, CImg<double> image1_pyramid, CImg<double> image2_pyramid)
{
	CImg<double> Blended_Laplace(image1_pyramid.height(), image1_pyramid.width(), image1_pyramid.depth(), image1_pyramid.spectrum());
	for(int c = 0; c < image1_pyramid.spectrum(); c++)
		for(int i = 0; i < image1_pyramid.height(); i++)
			for(int j = 0; j < image1_pyramid.width(); j++)
				Blended_Laplace.atXY(i, j, c) = (255 - mask_pyramid.atXY(i, j))*image1_pyramid.atXY(i, j, c) + mask_pyramid.atXY(i, j)*image2_pyramid.atXY(i, j, c);
	return Blended_Laplace;
}

CImg<double> Image_Addition(CImg<double> image1, CImg<double> image2)
{
	int rows1 = image1.height();
	int cols1 = image1.width();
	CImg<double> added_image(rows1, cols1, image1.depth(), image1.spectrum());
	for(int c = 0; c < 3; c++)
		for(int i = 0; i < rows1; i++)
			for(int j = 0; j < cols1; j++)
				added_image.atXY(i,j,c) = image1.atXY(i,j,c) + image2.atXY(i,j,c); 	
	return added_image;
}


int main(int argc, char **argv)
{
	string part = argv[1];
  try {
    

    //billboarrd1.png - 101,61 532,61 101,203 532,203
    if(part == "part1"){
      // Billboard
	      CImg<double> input_image("images/part1/lincoln.png");
	// Storing the transformation matrix
	CImg<double> transformation(3,3);
	transformation(0,0)=0.907;
	transformation(0,1)=0.258;
	transformation(0,2)=-182;
	transformation(1,0)=-0.153;
	transformation(1,1)=1.44;
	transformation(1,2)=58;
	transformation(2,0)=-0.000306;
	transformation(2,1)=0.000731;	
	transformation(2,2)=1;		
	
	//Creating a warped image by applying the transformation matrix
	CImg<double> transformed_image=apply_transformation(input_image,transformation);
	transformed_image.save("lincoln_warped.png");
	
    CImg<double> image1("images/part1/book1.jpg");
    CImg<double> image2("images/part1/book2.jpg");
	
	/////------------------------Homography solving------------------------------------------
	//Storing the four points for each images in their corresponding vectors
	double x_1=318,y_1=256,x1=141,y1=131;
	double x_2=534,y_2=372,x2=480,y2=159;
	double x_3=316,y_3=670,x3=493,y3=630;
	double x_4=73, y_4=473,x4=64, y4=601;	
	vector<Point> img1_points;
	img1_points.push_back({x_1,y_1});
	img1_points.push_back({x_2,y_2});
	img1_points.push_back({x_3,y_3});
	img1_points.push_back({x_4,y_4});
	vector<Point> img2_points;
	img2_points.push_back({x1,y1});
	img2_points.push_back({x2,y2});
	img2_points.push_back({x3,y3});
	img2_points.push_back({x4,y4});
	// Solving for homography given the above 4 matches
	CImg<double> homography_matrix= solve_homography(img1_points,img2_points);
	//Applying the homography on the second image
	CImg<double> warped_image2=apply_transformation(image2,homography_matrix);
	warped_image2.save("warped_image2.png");
	
	/////------------------------Billboard Creation------------------------------------------
	CImg<double> billboard1("images/part1/billboard1.jpg");
	vector<Point> billboard1_corners;
	billboard1_corners.push_back({101,61});
	billboard1_corners.push_back({532,61});
	billboard1_corners.push_back({101,203});
	billboard1_corners.push_back({532,203});
    CImg<double> billboard2("images/part1/billboard2.png");
	vector<Point> billboard2_corners;
	billboard2_corners.push_back({174,54});
	billboard2_corners.push_back({1107,260});
	billboard2_corners.push_back({148,623});
	billboard2_corners.push_back({1124,701});
	CImg<double> billboard3("images/part1/billboard3.jpg");
	vector<Point> billboard3_corners;
	billboard3_corners.push_back({618,287});
	billboard3_corners.push_back({1258,262});
	billboard3_corners.push_back({611,606});
	billboard3_corners.push_back({1260,601});
    char* poster_name = argv[2];
    CImg<double> poster(poster_name);
	CImg<double> billboard_poster_1=create_billboard(billboard1,billboard1_corners,poster);
	billboard_poster_1.save("synthetic_billboard1.png");
	CImg<double> billboard_poster_2=create_billboard(billboard2,billboard2_corners,poster);
	billboard_poster_2.save("synthetic_billboard2.png");
	CImg<double> billboard_poster_3=create_billboard(billboard3,billboard3_corners,poster);
	billboard_poster_3.save("synthetic_billboard3.png");
    }	
	  else if(part == "part2")
	  {
      char* image1 = argv[2];
	  char* image2 = argv[3];
	  char* image_mask_filename = argv[4];
	  // Image	  Blending 
	  // Resizing the images:
	  CImg<double> input_image1(image1);
	  CImg<double> image_1 = input_image1.get_resize(307,307);
	  CImg<double> input_image2(image2);
	  CImg<double> image_2 = input_image2.get_resize(307,307);
	  CImg<double> input_image_mask(image_mask_filename);
	  CImg<double> image_mask = input_image_mask.get_resize(307,307);
	  
	  // Defining the Gaussian Kernel:
	  CImg<double> Gaussian_kernel(5,5);
	  Gaussian_kernel(0,0) = float(1)/float(256); Gaussian_kernel(0,1) = float(4)/float(256); Gaussian_kernel(0,2) = float(6)/float(256); Gaussian_kernel(0,3) = float(4)/float(256); Gaussian_kernel(0,4) = float(1)/float(256);
	  
	  Gaussian_kernel(1,0) = float(4)/float(256); Gaussian_kernel(1,1) = float(16)/float(256); Gaussian_kernel(1,2) = float(24)/float(256); Gaussian_kernel(1,3) = float(16)/float(256); Gaussian_kernel(1,4) = float(4)/float(256);
	  
	  Gaussian_kernel(2,0) = float(6)/float(256); Gaussian_kernel(2,1) = float(24)/float(256); Gaussian_kernel(2,2) = float(36)/float(256); Gaussian_kernel(2,3) = float(24)/float(256); Gaussian_kernel(2,4) = float(6)/float(256);
	  
	  Gaussian_kernel(3,0) = float(4)/float(256); Gaussian_kernel(3,1) = float(16)/float(256); Gaussian_kernel(3,2) = float(24)/float(256); Gaussian_kernel(3,3) = float(16)/float(256); Gaussian_kernel(3,4) = float(4)/float(256);
	  
	  Gaussian_kernel(4,0) = float(1)/float(256); Gaussian_kernel(4,1) = float(4)/float(256); Gaussian_kernel(4,2) = float(6)/float(256); Gaussian_kernel(4,3) = float(4)/float(256); Gaussian_kernel(4,4) = float(1)/float(256);
			  
	  CImg<double> Gmag4_kernel(5,5);
	  Gmag4_kernel(0,0) = float(1)/float(64); Gmag4_kernel(0,1) = float(4)/float(64); Gmag4_kernel(0,2) = float(6)/float(64); 
	  Gmag4_kernel(0,3) = float(4)/float(64); Gmag4_kernel(0,4) = float(1)/float(64);
	  
	  Gmag4_kernel(1,0) = float(4)/float(64); Gmag4_kernel(1,1) = float(16)/float(64); Gmag4_kernel(1,2) = float(24)/float(64); 
	  Gmag4_kernel(1,3) = float(16)/float(64); Gmag4_kernel(1,4) = float(4)/float(64);
	  
	  Gmag4_kernel(2,0) = float(6)/float(64); Gmag4_kernel(2,1) = float(24)/float(64); Gmag4_kernel(2,2) = float(36)/float(64); 
	  Gmag4_kernel(2,3) = float(24)/float(64); Gmag4_kernel(2,4) = float(6)/float(64);
	  
	  Gmag4_kernel(3,0) = float(4)/float(64); Gmag4_kernel(3,1) = float(16)/float(64); Gmag4_kernel(3,2) = float(24)/float(64); 
	  Gmag4_kernel(3,3) = float(16)/float(64); Gmag4_kernel(3,4) = float(4)/float(64);
	  
	  Gmag4_kernel(4,0) = float(1)/float(64); Gmag4_kernel(4,1) = float(4)/float(64); Gmag4_kernel(4,2) = float(6)/float(64); 
	  Gmag4_kernel(4,3) = float(4)/float(64); Gmag4_kernel(4,4) = float(1)/float(64);
	  
	  //Gaussian Pyramid for image 1:
	  CImg<double> g0_img1 = image_1;	  
	  CImg<double> g1_img1 = Gaussian_Pyramid(g0_img1, Gaussian_kernel);	  
	  CImg<double> g2_img1 = Gaussian_Pyramid(g1_img1, Gaussian_kernel);	  
	  CImg<double> g3_img1 = Gaussian_Pyramid(g2_img1, Gaussian_kernel);	  
	  CImg<double> g4_img1 = Gaussian_Pyramid(g3_img1, Gaussian_kernel);	  
	  CImg<double> g5_img1 = Gaussian_Pyramid(g4_img1, Gaussian_kernel);	  
	  
	  //Laplacian of Gaussian for image 1:
	  CImg<double> lg1_img1 = g1_img1.get_resize(g0_img1.width(), g0_img1.height());	  
	  CImg<double> lg2_img1 = g2_img1.get_resize(g1_img1.width(), g1_img1.height());		  
	  CImg<double> lg3_img1 = g3_img1.get_resize(g2_img1.width(), g2_img1.height());	  
	  CImg<double> lg4_img1 = g4_img1.get_resize(g3_img1.width(), g3_img1.height());	  
	  CImg<double> lg5_img1 = g5_img1.get_resize(g4_img1.width(), g4_img1.height());	  
	  
	  //Laplacian Pyramid for image 1:
	  CImg<double> l0_img1 = Difference(g0_img1, lg1_img1); 
	  CImg<double> l1_img1 = Difference(g1_img1, lg2_img1); 
	  CImg<double> l2_img1 = Difference(g2_img1, lg3_img1); 
	  CImg<double> l3_img1 = Difference(g3_img1, lg4_img1); 
	  CImg<double> l4_img1 = Difference(g4_img1, lg5_img1); 
	  CImg<double> l5_img1 = g5_img1; 
	  
	  //Gaussian Pyramid for image 2:
	  CImg<double> g0_img2 = image_2;
	  CImg<double> g1_img2 = Gaussian_Pyramid(g0_img2, Gaussian_kernel);
	  CImg<double> g2_img2 = Gaussian_Pyramid(g1_img2, Gaussian_kernel);
	  CImg<double> g3_img2 = Gaussian_Pyramid(g2_img2, Gaussian_kernel);
	  CImg<double> g4_img2 = Gaussian_Pyramid(g3_img2, Gaussian_kernel);
	  CImg<double> g5_img2 = Gaussian_Pyramid(g4_img2, Gaussian_kernel);
	  
	  //Laplacian of Gaussian for image 2:
	  CImg<double> lg1_img2 = g1_img2.get_resize(g0_img2.width(), g0_img2.height());
	  CImg<double> lg2_img2 = g2_img2.get_resize(g1_img2.width(), g1_img2.height());
	  CImg<double> lg3_img2 = g3_img2.get_resize(g2_img2.width(), g2_img2.height());
	  CImg<double> lg4_img2 = g4_img2.get_resize(g3_img2.width(), g3_img2.height());
	  CImg<double> lg5_img2 = g5_img2.get_resize(g4_img2.width(), g4_img2.height());
	  
	  //Laplacian Pyramid for image 2:
	  CImg<double> l0_img2 = Difference(g0_img2, lg1_img2);
	  CImg<double> l1_img2 = Difference(g1_img2, lg2_img2);
	  CImg<double> l2_img2 = Difference(g2_img2, lg3_img2);
	  CImg<double> l3_img2 = Difference(g3_img2, lg4_img2);
	  CImg<double> l4_img2 = Difference(g4_img2, lg5_img2);
	  CImg<double> l5_img2 = g5_img2;
	  
	  //Gaussian Pyramid for mask:
	  CImg<double> g0_mask = image_mask;
	  CImg<double> g1_mask = Gaussian_Pyramid(g0_mask, Gaussian_kernel);
	  CImg<double> g2_mask = Gaussian_Pyramid(g1_mask, Gaussian_kernel);
	  CImg<double> g3_mask = Gaussian_Pyramid(g2_mask, Gaussian_kernel);
	  CImg<double> g4_mask = Gaussian_Pyramid(g3_mask, Gaussian_kernel);
	  CImg<double> g5_mask = Gaussian_Pyramid(g4_mask, Gaussian_kernel);
	  
	  // Blended Laplacian Pyramids
	  CImg<double> bl0 = Blended_Laplacian(g0_mask, l0_img1, l0_img2);
	  CImg<double> bl1 = Blended_Laplacian(g1_mask, l1_img1, l1_img2);
	  CImg<double> bl2 = Blended_Laplacian(g2_mask, l2_img1, l2_img2);
	  CImg<double> bl3 = Blended_Laplacian(g3_mask, l3_img1, l3_img2);
	  CImg<double> bl4 = Blended_Laplacian(g4_mask, l4_img1, l4_img2);
	  CImg<double> bl5 = Blended_Laplacian(g5_mask, l5_img1, l5_img2);
	  
	  // Image Reconstruction:
	  CImg<double> f5_1 = bl5.get_resize(19, 19);
	  CImg<double> f5 = f5_1.get_convolve(Gaussian_kernel);
	  
	  CImg<double> f4_1 = Image_Addition(f5,bl4);
	  CImg<double> f4_2 = f4_1.get_resize(38, 38);
	  CImg<double> f4 = f4_2.get_convolve(Gaussian_kernel);
	  
	  CImg<double> f3_1 = Image_Addition(f4,bl3);
	  CImg<double> f3_2 = f3_1.get_resize(76, 76);
	  CImg<double> f3 = f3_2.get_convolve(Gaussian_kernel);
	  
	  CImg<double> f2_1 = Image_Addition(f3,bl2);
	  CImg<double> f2_2 = f2_1.get_resize(153, 153);
	  CImg<double> f2 = f2_2.get_convolve(Gaussian_kernel);
	  
	  CImg<double> f1_1 = Image_Addition(f2,bl1);
	  CImg<double> f1_2 = f1_1.get_resize(307, 307);
	  CImg<double> f1 = f1_2.get_convolve(Gaussian_kernel);
	  
	  CImg<double> f0 = Image_Addition(f1,bl0);
	  f0.save("Blended_Output.png");
	  }

    else if(part == "part3"){
	char* image1_filename = argv[2];
	char* image2_filename = argv[3];
	CImg<double> part3_image1(image1_filename);
    CImg<double> part3_image2(image2_filename);
    vector<PointPair> matches = find_matches_in_images(part3_image1,part3_image2,0.9);
	CImg<double> good_match_image=draw_lines_sift_matches(part3_image1,part3_image2,matches);
	good_match_image.save("good_match_image.png");
	struct homogdist  ransac_output=ransac(matches,3.0);
	CImg<double> ransac_match_image=draw_lines_sift_matches(part3_image1,part3_image2,ransac_output.good_matches);
	ransac_match_image.save("ransac_match_image.png");
	// RANSAC
    }
    else if(part == "part4"){

	char* image1_filename = argv[2];
	  char* image2_filename = argv[3];
	  char* image3_filename = argv[4];
		CImg<double> part4_image1(image1_filename);
		CImg<double> part4_image2(image2_filename);
		CImg<double> part4_image3(image3_filename);
     	vector<trans_images_struct> transform_images_out=transform_images(part4_image1,part4_image2,part4_image3);
		vector<const char*> names;
		names.push_back("warped_image_1.png");
		names.push_back("warped_image_3.png");
		names.push_back("panorama.png");
		//for(int i=0;i<transform_images_out.size();i++)
		transform_images_out[2].image.save(names[2]);

    }
    
    
    // feel free to add more conditions for other parts (e.g. more specific)
    //  parts, for debugging, etc.
  }
  catch(const string &err) {
    cerr << "Error: " << err << endl;
  }
}
