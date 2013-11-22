#include "itkImage.h"
#include "itkImageFileReader.h"
#include "itkImageFileWriter.h"
#include "itkParticle.h"

//------------------------------------------------------------------------------
//
// ./ParticleFilter FixedImage MovingImage TransformedImage DisplacementField InputDisplacement
//
// FixedImage                     -> Fixed image file name (Input)
// MovingImage                    -> Moving image file name (Input)
// TransformedImage         			<- Name of file to save transformed image (Output)
// DisplacementField              <- Name of Displacement Field (Output Optional)
//
//------------------------------------------------------------------------------

int main( int argc, char **argv ) {
        //Types definition
        const int        Dimensions = 3;
        typedef itk::Image< short, Dimensions >													ImageType;
        typedef itk::Image< float, Dimensions >                         FloatImageType;
        typedef itk::Vector<float, Dimensions>                          VectorPixelType;
        typedef itk::Image<VectorPixelType,Dimensions>                  VectorImageType;
        typedef itk::ImageFileReader< ImageType >                       ReaderType;
        typedef itk::ImageFileReader< VectorImageType >                 ReaderVectorType;
        typedef itk::ImageFileWriter< ImageType >                       WriterType;
        typedef itk::Particle<ImageType,ImageType >                  		ParticleType;
        typedef itk::ImageFileWriter< VectorImageType >                 WriterVectorType;

   /*     //Create Objects
        ReaderType::Pointer                                 fixedImage = ReaderType::New();
        ReaderType::Pointer                                 movingImage = ReaderType::New();
        ReaderVectorType::Pointer                vectorImage = ReaderVectorType::New();
        WriterType::Pointer                                 writer = WriterType::New();
        WriterVectorType::Pointer                vectorWriter = WriterVectorType::New();
        OpticalFlowType::Pointer                opticalFlow = OpticalFlowType::New();

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

        //Set variables to optical flow
        opticalFlow->SetFixedImage( fixedImage->GetOutput() );
        opticalFlow->SetMovingImage( movingImage->GetOutput() );

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
*/
        return 0;
}


