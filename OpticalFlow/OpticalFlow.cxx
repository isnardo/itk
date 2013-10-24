#include "itkImage.h"
#include "itkImageFileReader.h"
#include "itkImageFileWriter.h"
#include "itkOpticalFlow.h"

//------------------------------------------------------------------------------
//
// ./OpticalFlow FixedImage MovingImage TransformedImage VectorField InputVectorField
// 
// FixedImage 				-> Fixed image file name (Input)
// MovingImage 				-> Moving image file name (Input)
// TransformedImage 	<- Name of file to save transformed image (Output)
// VectorField 				<- Name of Vector Field (Output Optional)
// VectorField 				-> Name of Input Vector Field (Input Optional)
//
//------------------------------------------------------------------------------


int main( int argc, char **argv ) {
	const int	Dimensions = 2;

	typedef itk::Image< unsigned char, Dimensions >				ImageType;
	typedef itk::Image< float, Dimensions >    						FloatImageType;
	typedef itk::Vector<float, Dimensions>								VectorPixelType;
	typedef itk::Image<VectorPixelType,Dimensions>				VectorImageType;
	typedef itk::ImageFileReader< ImageType >							ReaderType;
	typedef itk::ImageFileReader< VectorImageType >				ReaderVectorType;
	typedef itk::ImageFileWriter< ImageType >							WriterType;
	typedef itk::OpticalFlow<ImageType,ImageType >				OpticalFlowType;
	typedef itk::ImageFileWriter< VectorImageType >	 			WriterVectorType;

	ReaderType::Pointer 				fixedImage = ReaderType::New();
	ReaderType::Pointer 				movingImage = ReaderType::New();
	ReaderVectorType::Pointer		vectorImage = ReaderVectorType::New();
	WriterType::Pointer 				writer =  WriterType::New();
	WriterVectorType::Pointer		vectorWriter = WriterVectorType::New();
	OpticalFlowType::Pointer		Filter = OpticalFlowType::New();


	return 0;
}
