#include "itkImage.h"
#include "itkImageFileReader.h"
#include "itkImageFileWriter.h"
#include "itkMultiScaleOpticalFlow.h"

//------------------------------------------------------------------------------
//
// ./OpticalFlow FixedImage MovingImage TransformedImage DisplacemntField
//
// FixedImage               -> Fixed image file name (Input)
// MovingImage              -> Moving image file name (Input)
// TransformedImage         <- Name of file to save transformed image (Output)
// DisplacemntField         <- Name of Vector Field (Output Optional)
//
//------------------------------------------------------------------------------

int main( int argc, char **argv ) {
        //Types definition
        const int        Dimensions = 3;
        typedef itk::Image< short, Dimensions >                               ImageType;
        typedef itk::Image< float, Dimensions >                               FloatImageType;
        typedef itk::Vector<float, Dimensions>                                VectorPixelType;
        typedef itk::Image<VectorPixelType,Dimensions>                        VectorImageType;
        typedef itk::ImageFileReader< ImageType >                             ReaderType;
        typedef itk::ImageFileWriter< ImageType >                             WriterType;
        typedef itk::MultiScaleOpticalFlow<ImageType,ImageType >              OpticalFlowType;
        typedef itk::ImageFileWriter< VectorImageType >                       WriterVectorType;

        //Create Objects
        ReaderType::Pointer                           fixedImage = ReaderType::New();
        ReaderType::Pointer                           movingImage = ReaderType::New();
        WriterType::Pointer                           writer = WriterType::New();
        WriterVectorType::Pointer                			vectorWriter = WriterVectorType::New();
        OpticalFlowType::Pointer                			opticalFlow = OpticalFlowType::New();

        //Load fixed image
        fixedImage->SetFileName( argv[1] );
        try
        {
                fixedImage->Update();
        }
        catch( itk::ExceptionObject & excep )
        {
                std::cerr << "Exception catched !" << std::endl;
                std::cerr << excep << std::endl;
        }

        //Load moving image
        movingImage->SetFileName( argv[2] );
        try
        {
                movingImage->Update();
        }
        catch( itk::ExceptionObject & excep )
        {
                std::cerr << "Exception catched !" << std::endl;
                std::cerr << excep << std::endl;
        }

				//Itrations for each scale
				int *Iterations = new int [4];

				Iterations[0] = 20; //Default 20
				Iterations[1] = 20; //Default 20
				Iterations[2] = 20;	//Default 20
				Iterations[3] = 2;	//Default 20

        //Set variables to optical flow
        opticalFlow->SetFixedImage( fixedImage->GetOutput() );
        opticalFlow->SetMovingImage( movingImage->GetOutput() );
        
        opticalFlow->SetLambdaXY( 1 );        //Default 1
        opticalFlow->SetLambdaZ( 1 );         //Default 1
        opticalFlow->SetLambdaW( 1 );         //Default 1
        opticalFlow->SetThreshold( 0.1 );     //Default 1.0

				opticalFlow->SetIterationsOfEachScale( Iterations );

				opticalFlow->ShowIterationsOn();	//Show Iterations at each scale (Default Off)
				opticalFlow->ShowLevelsOn();			//Show error at each level			(Default Off)

				//Update optical flow
        try
        {
                opticalFlow->Update();
        }
        catch( itk::ExceptionObject & excep )
        {
                std::cerr << "Exception catched !" << std::endl;
                std::cerr << excep << std::endl;
        }

				std::cout<<"LLEGO"<<std::endl;
       //save image transformed
        writer->SetFileName( argv[3] );
        writer->SetInput( opticalFlow->GetTransformedImage() );
        try
        {
                writer->Update();
        }
        catch( itk::ExceptionObject & excep )
        {
                std::cerr << "Exception catched !" << std::endl;
                std::cerr << excep << std::endl;
        }

        //save displacements vector field
        if( argc > 4)
        {
        	vectorWriter->SetFileName( argv[4] );
        	vectorWriter->SetInput( opticalFlow->GetDisplacementField() );
        	try
        	{
        		vectorWriter->Update();
          }
          catch( itk::ExceptionObject & excep )
          {
          	std::cerr << "Exception catched !" << std::endl;
            std::cerr << excep << std::endl;
          }
        }

        return 0;
}


