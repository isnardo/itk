#ifndef _itkOpticalFlow_h
#define _itkOpticalFlow_h

#include "itkImage.h"
#include "itkVector.h"
#include "itkImageToImageFilter.h"
#include "itkImageDuplicator.h"
#include "itkCastImageFilter.h"
#include "itkAddImageFilter.h"
#include "itkSubtractImageFilter.h"
#include "itkConvolutionImageFilter.h"
#include "itkLinearInterpolateImageFunction.h"
#include "itkMinimumMaximumImageCalculator.h"
#include "itkWarpImageFilter.h"

namespace itk
{
	//Template Class Optical Flow
	template< class TInputImage, class TOutputImage>
	class ITK_EXPORT OpticalFlow : public ImageToImageFilter< TInputImage, TOutputImage >
	{
		public:
			//** Standard class typedefs **//
			typedef OpticalFlow      																Self;
			typedef ImageToImageFilter< TInputImage,TOutputImage >	Superclass;
			typedef SmartPointer< Self >        										Pointer;
			typedef SmartPointer< const Self >        							ConstPointer;

			itkNewMacro( Self );
			itkTypeMacro( OpticalFlow,ImageToImageFilter );

			//Types definition
			typedef TInputImage 																		InputImageType;
			typedef TOutputImage 																		OutputImageType;

			//Const Macro to Define the Images Dimensions
			itkStaticConstMacro(ImageDimension, unsigned int,InputImageType::ImageDimension);
			//Type to manage an internal image type
			typedef itk::Image< float,itkGetStaticConstMacro(ImageDimension)>		
																															InternalImageType;
			//Input Caster
 			typedef itk::CastImageFilter<InputImageType,InternalImageType >
																															InputCasterType;
			//Output Caster
 			typedef itk::CastImageFilter<InternalImageType,OutputImageType>
																															OutputCasterType;
			//Interpolator
			typedef InterpolateImageFunction<InternalImageType,double> 	
																															InterpolatorType;
			//Default Interpolator Type (Linear Interpolator)
			typedef LinearInterpolateImageFunction<InternalImageType,double>
																															DefaultInterpolatorType;
			//Vector pixel type for the displacements vector field
			typedef itk::Vector< float,itkGetStaticConstMacro(ImageDimension)>	
																															VectorPixelType;
			//Vector Image Type for the displacements vector field
			typedef itk::Image< VectorPixelType,itkGetStaticConstMacro(ImageDimension)>								
																															VectorImageType;
			//Convolution filter to compute image derivates
			typedef itk::ConvolutionImageFilter<InternalImageType,InternalImageType,InternalImageType>					
																															ConvolutionFilterType;
			//Add filter (+)
			typedef itk::AddImageFilter<VectorImageType,VectorImageType>	
																															AddFilterType;
			//Substraction filter (-)
			typedef itk::SubtractImageFilter<InternalImageType,InternalImageType,InternalImageType>		
																															SubFilterType;
			//Duplicator Filter to clone images
			typedef itk::ImageDuplicator< VectorImageType > 				DuplicatorType;
			//Warp image filter to apply displacements vector field
			typedef typename itk::WarpImageFilter<InternalImageType,InternalImageType,VectorImageType>
																															WarpFilterType;
			//Min-Max filter
			typedef itk::MinimumMaximumImageCalculator<InternalImageType>  
																															ImageMinMaxFilterType;

			//Macro to Set other type of interpolator
			itkSetObjectMacro( Interpolator,InterpolatorType );
			//Iterations number for the optical flow 
			itkSetMacro( m_Iterations,int );
			itkGetMacro( m_Iterations,int );
			//Iterations number for Gauss-Seidel method to solve the linear sistem of equations
			itkSetMacro( m_GaussSeidelIterations,int );
			itkGetMacro( m_GaussSeidelIterations,int );
			//Stop Threshold according to intensities differences
			itkSetMacro( m_Threshold,float );
			itkGetMacro( m_Threshold,float );
			//Lambda XY (Weigth in directions X and Y for the regularization term in the energy function)
			itkSetMacro( m_LambdaXY,float );
			itkGetMacro( m_LambdaXY,float );
			//Lambda Z (Weigth in directions Z for the regularization term in the energy function)
			itkSetMacro( m_LambdaZ,float );
			itkGetMacro( m_LambdaZ,float );
			//Lambda W (Weigth for all displacements)
			itkSetMacro( m_LambdaW,float );
			itkGetMacro( m_LambdaW,float );


			//Types image elements
			typedef typename InternalImageType::PixelType						PixelType;
			typedef typename InternalImageType::IndexType						IndexType;
			typedef typename InternalImageType::SizeType						SizeType;
			typedef typename InternalImageType::RegionType					RegionType;
			typedef typename VectorImageType::IndexType							VectorIndexType;

			
			//** Public functions definition **//
			virtual void SetFixedImage ( const  InputImageType  *fixedImage  );
			virtual void SetMovingImage( const  InputImageType  *movingImage );
			virtual void SetInternalFixedImage( InternalImageType  *fixedImage )
										{
											FixedImage = fixedImage;
											FixedImageSize = FixedImage->GetLargestPossibleRegion().GetSize();
										};
			virtual void SetInternalMovingImage( InternalImageType *movingImage )
										{
											MovingImage = movingImage;
											MovingImageSize = MovingImage->GetLargestPossibleRegion().GetSize();
										};
			virtual void SetDisplacementField(typename VectorImageType::Pointer DisplacementField_t)
										{
											DisplacementField = DisplacementField_t;
											DisplacementFieldFlag = true;
										};
			virtual typename DisplacementImageType::Pointer GetDisplacementField()
												{
													return DisplacementField;
												};
			virtual typename OutputImageType::Pointer GetTransformedImage();

		protected:

			OpticalFlow();
			~OpticalFlow(){};

			virtual void GenerateData();

			//Transform Internal Image Type using a Given a Displacement Field
			virtual typename InternalImageType::Pointer ApplyDisplacementField(typename InternalImageType::Pointer I,
																							 typename VectorImageType::Pointer V)
			{

				typename ImageMinMaxFilterType::Pointer		MinMax  = ImageMinMaxFilterType::New();
				typename WarpFilterType::Pointer 					WarpFilter = WarpFilterType::New();

				//Find the Minimum Value
				MinMax->SetImage( I );
				MinMax->Compute();
				float minimum = MinMax->GetMinimum();	

				//Set Warp Filter Parameters
				WarpFilter->SetDisplacementField( V );
				WarpFilter->SetInterpolator( m_Interpolator );
				WarpFilter->SetInput( I );
				WarpFilter->SetOutputOrigin( I->GetOrigin() );
				WarpFilter->SetOutputSpacing( I->GetSpacing() );
				WarpFilter->SetOutputDirection( I->GetDirection() );
				WarpFilter->SetEdgePaddingValue( minimum );
				//Update Warp Filter
				try
				{
					WarpFilter->Update();
				}
				catch( itk::ExceptionObject & excep )
				{
					std::cerr << "Exception catched !" << std::endl;
					std::cerr << excep << std::endl;
				}
			
				//Return the deformed Image
				return WarpFilter->GetOutput();
			};

		private:
			//Private Variables
			float																				m_LambdaXY;
			float																				m_LambdaZ;
			float																				m_LambdaW;

			int 																				m_Iterations;
			int 																				m_GaussSeidelIterations;
			float																				m_Threshold;

			bool																				DisplacementFieldFlag;
			float																				IntensitiesRange;

			SizeType																		FixedImageSize;
			SizeType																		MovingImageSize;

			typename InternalImageType::Pointer					MovingImage;
			typename InternalImageType::Pointer					FixedImage;
			typename InputImageType::ConstPointer				InputFixedImage;
			typename InputImageType::ConstPointer				InputMovingImage;
			typename InterpolatorType::Pointer					m_Interpolator;
			typename InternalImageType::Pointer					DerivativeX;
			typename InternalImageType::Pointer					DerivativeY;
			typename InternalImageType::Pointer					DerivativeZ;
			typename InternalImageType::Pointer					aux_DerivativeX;
			typename InternalImageType::Pointer					aux_DerivativeY;
			typename InternalImageType::Pointer					aux_DerivativeZ;
			typename InternalImageType::Pointer					SubImage;
			typename VectorImageType::Pointer						DisplacementField;
			typename VectorImageType::Pointer						aux_DeformationField;

			void 	ComputeGradient();
			void 	SolveGaussSeidel3D();
			void 	SolveGaussSeidel2D();
			void 	ComputeOpticalFlow();
			void 	AllocateDisplacementField();
			float ComputeIntensitiesRange();
			float ComputePercentageImageDifferences();

	};//Class Optical Flow

	
} //namespace ITK

#ifndef ITK_MANUAL_INSTANTIATION
#include "itkOpticalFlow.txx"
#endif

#endif