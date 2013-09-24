#include "itkImage.h"
#include "itkImageFileReader.h"
#include "itkImageFileWriter.h"
#include "itkParticle.h"
#include "itkParticleFilter.h"
#include "LinuxTimer.h"
#include "itkBSplineInterpolateImageFunction.h"
#include "itkLinearInterpolateImageFunction.h"
#include "itkNearestNeighborInterpolateImageFunction.h"
#include "itkWindowedSincInterpolateImageFunction.h"
//------------------------------------------------------------------------------
//
// ./PF FixedImage MovingImage TransformedImage DisplacementField
// 
// FixedImage 				-> Fixed image file name (Input)
// MovingImage 			-> Moving image file name (Input)
// TransformedImage 		<- Name of file to save transformed image (Output)
// DisplacementField 	<- Name of Vector Field (Output Optional)
//
//------------------------------------------------------------------------------

int main( int argc, char **argv ) {
	const int	Dimensions = 3;
	typedef itk::Image< unsigned char, Dimensions >					ImageType;
//	typedef itk::Image< unsigned short, Dimensions >				ImageType;
	typedef itk::Image< float, Dimensions >    						FloatImageType;
	typedef itk::Vector<float, Dimensions>								VectorPixelType;
	typedef itk::Image<VectorPixelType,Dimensions>					VectorImageType;
	typedef itk::Image<float,Dimensions >    							FloatImageType;
	typedef itk::ImageFileReader< ImageType >							ReaderType;
	typedef itk::ImageFileWriter< ImageType >			 				WriterType;
	typedef itk::ImageFileWriter< VectorImageType >	 				WriterVectorType;
	typedef itk::ParticleFilter<ImageType,ImageType >				ParticleFilterType;
	typedef ParticleFilterType::Particle::VectorType				VectorType;

	Linux::Timer t;	//Timer

	ReaderType::Pointer 				fixedImage = ReaderType::New();
	ReaderType::Pointer 				movingImage = ReaderType::New();
	ParticleFilterType::Pointer	Filter =  ParticleFilterType::New();
	WriterType::Pointer 				writer =  WriterType::New();
	WriterVectorType::Pointer		vectorWriter = WriterVectorType::New();
	VectorType							Variances;
	VectorType							Parameters;

	//Load Images
	fixedImage->SetFileName( argv[1] );
	movingImage->SetFileName( argv[2] );
	try{
		fixedImage->Update();
		movingImage->Update();
	}
	catch( itk::ExceptionObject & excep )
	{
		std::cerr << "Exception catched !" << std::endl;
		std::cerr << excep << std::endl;
	}

	// Define number of parameters
	// 
	int NumberOfParameters = 9;

	Variances.SetSize( NumberOfParameters );
	Parameters.SetSize( NumberOfParameters );

/*	Parameters[0] = 0.0 ;
	Parameters[1] = 0.0 ;
	Parameters[2] = 0.0 ;
	Parameters[3] = 0.0 ;
	Parameters[4] = -40.0 ;
	Parameters[5] = 0.0 ;
	Parameters[6] = 1.5 ;
	Parameters[7] = 1.5 ;
	Parameters[8] = 1.0 ;

	Variances[0] = 0.1;	//Angle x in Degrees
	Variances[1] = 0.1;	//Angle y in Degrees
	Variances[2] = 1.0;	//Angle z in Degrees
	Variances[3] = 15.0;	//X Shift
	Variances[4] = 15.0;	//Y Shift
	Variances[5] = 1.0;	//Z Shift
	Variances[6] = 2.5e-3;	//X Scale
	Variances[7] = 2.5e-3;	//Y Scale  
	Variances[8] = 2.5e-3;	//Z Scale
*/


	ImageType::IndexType		Step;
	Step[0] = 8;
	Step[1] = 8;
	Step[2] = 2;


	Filter->SetFixedImage( fixedImage->GetOutput() );
	Filter->SetMovingImage( movingImage->GetOutput() );
//	Filter->SetNumberOfParameters( NumberOfParameters );
//	Filter->SetIntialParameters( Parameters );
	Filter->SetVarianceOfPerturbations( Variances );
	Filter->SetNumberOfParticles( 400 );
	Filter->SetMaxIterations( 400 );
	Filter->SetSamplesStepSize( Step );
	Filter->SetNumberOfBins( 64 );
//	Filter->PhysicalCoordinatesOn(); //Use the physical coordinates (default use voxel coordinates)
//	Filter->SetInterpolator( Interpolator );


	t.start();
	try{
		Filter->Update();
	}
	catch( itk::ExceptionObject & excep )
	{
		std::cerr << "Exception catched !" << std::endl;
		std::cerr << excep << std::endl;
	}
	t.stop();

	std::cout<<Filter->GetParameters()<<std::endl;	
	std::cout<<t.duration()/60000.0<<" min"<<std::endl;

	writer->SetFileName( argv[3] );
	writer->SetInput( movingImage->GetOutput() );
	writer->SetInput( Filter->GetTransformedImage() );
	writer->Update();

	if( argc > 4)
	{
		vectorWriter->SetFileName( argv[4] );
		vectorWriter->SetInput( Filter->GetDisplacementField() );
		vectorWriter->Update();
	}


	return 0;
}
