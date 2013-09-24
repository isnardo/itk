#ifndef _itkParticle_h
#define _itkParticle_h

#include "itkImage.h"
#include "itkImageToImageFilter.h"
#include "itkVariableLengthVector.h"
#include "itkVector.h"
#include "itkMinimumMaximumImageCalculator.h"
#include "itkImageRegionIterator.h"
#include "itkImageRegionIteratorWithIndex.h"
#include "itkCastImageFilter.h"
#include "itkLinearInterpolateImageFunction.h"
#include "random.h"
#include "itkImageDuplicator.h"
#include "itkContinuousIndex.h"
#include "itkWarpImageFilter.h"

namespace itk
{
	//Particle Class
	template< class TInputImage, class TOutputImage>
	class ITK_EXPORT Particle : public ImageToImageFilter< TInputImage, TOutputImage >
	{
		public:
			// Standard class typedefs
			typedef Particle             												Self;
			typedef ImageToImageFilter< TInputImage,TOutputImage >			Superclass;
			typedef SmartPointer< Self >        									Pointer;
			typedef SmartPointer< const Self >        							ConstPointer;

			itkNewMacro(Self);
			itkTypeMacro(Particle, ImageToImageFilter);

			//Types definition
			typedef TInputImage 															InputImageType;
			typedef TOutputImage 														OutputImageType;

			//Const Macro to Define the Images Dimensions
			itkStaticConstMacro(ImageDimension, unsigned int,InputImageType::ImageDimension);

			typedef itk::Image< float,itkGetStaticConstMacro(ImageDimension)>		
																								InternalImageType;
			typedef itk::Vector< float,itkGetStaticConstMacro(ImageDimension)>	
																								VectorPixelType;
			typedef itk::Image< VectorPixelType,
										itkGetStaticConstMacro(ImageDimension)>								
																								VectorImageType;
			typedef itk::Matrix< float,4,4 > 										MatrixType;
			typedef itk::VariableLengthVector< float > 							VectorType;
			typedef itk::MinimumMaximumImageCalculator<InternalImageType>  ImageMinMaxFilterType;
 			typedef itk::CastImageFilter<InputImageType,InternalImageType >  																									InputCasterType;
 			typedef itk::CastImageFilter<InternalImageType,OutputImageType >  																								OutputCasterType;

			typedef double 																CoordRepType;
			typedef InterpolateImageFunction<InternalImageType,CoordRepType> 	
																								InterpolatorType;
			typedef LinearInterpolateImageFunction<InternalImageType,CoordRepType >
																								DefaultInterpolatorType;
			typedef itk::ImageDuplicator< InternalImageType > 					InternalDuplicatorType;

			typedef typename InternalImageType::PixelType						PixelType;
			typedef typename InternalImageType::IndexType						IndexType;
			typedef typename InternalImageType::SizeType							SizeType;
			typedef typename InternalImageType::PointType						PointType;
			typedef typename VectorImageType::IndexType							VectorIndexType;
	
			typedef typename itk::ContinuousIndex<double,itkGetStaticConstMacro(ImageDimension)> 
																								ContinuousIndex;
			typedef typename itk::MinimumMaximumImageCalculator<InternalImageType>
																								MinMaxFilterType;
			typedef typename itk::WarpImageFilter<InternalImageType,InternalImageType,VectorImageType>
																								WarpFilterType;

			itkSetObjectMacro( Interpolator, InterpolatorType );
			itkSetMacro( PhysicalCoordinates, bool);
			itkBooleanMacro( PhysicalCoordinates );


			//Public Functions
			virtual void SetFixedImage (const  InputImageType  *fixedImage  );
			virtual void SetMovingImage(const  InputImageType *movingImage );
			virtual void SetInternalFixedImage (const  InternalImageType  *fixedImage )
										{
											FixedImage = fixedImage;
											FixedImageSize = FixedImage->GetLargestPossibleRegion().GetSize();
											FixedSpacing =  FixedImage->GetSpacing();
										}
			virtual void SetInternalMovingImage(const  InternalImageType *movingImage )
										{
											MovingImage = movingImage;
											MovingImageSize = MovingImage->GetLargestPossibleRegion().GetSize();
											MovingSpacing = MovingImage->GetSpacing();
										}
			virtual void SetNumberOfParameters(int NumberOfParameters_t);
			virtual void SetNoiseStandardDeviation(float StdDev);
			virtual void SetSamplesStepSize(IndexType step){SamplesStepSize = step;};
			virtual void SetParameters(VectorType Vector);
			virtual void SetVarianceOfPerturbations(VectorType Vector);

			virtual void ReducePerturbations(float facDeg);
			void UpsetParticle();
			void UpsetParticle(VectorType Perturbations);
			void GenerateDisplacementField();
			float ComputeParticleWeight();
			void ScaleAddVector(VectorType Vec, float fac);

			virtual void SetNumberOfBins(int NumberOfBins_t)
								{
									NumberOfBins = NumberOfBins_t;
								};
			void SetParticleWeight(float Weight_t)
								{
									Weight = Weight_t;
								};
			float GetParticleWeight()
								{
									return Weight;
								};
			VectorType GetParameters()
							{
								return ParticleParameters;
							};
			virtual void SetPerturbations(VectorType Perturbations_t)
								{
									ParticlePerturbations = Perturbations_t;
								};
			virtual void FillParameters(float val)
								{
									ParticleParameters.Fill(val);
								};

			virtual typename OutputImageType::Pointer GetTransformedImage();
			virtual typename OutputImageType::Pointer GetTransformedImage(
																	typename VectorImageType::Pointer Vec_t)
																	{
																		DisplacementField = Vec_t;
																		ban_send_vec = true;
																		return GetTransformedImage();
																	};

			virtual typename VectorImageType::Pointer GetDisplacementField(){
																		GenerateDisplacementField();
																		return DisplacementField;
																	};

		protected:
			Particle();//Builder
			~Particle(){};//Destroy

			//Main Function of the Class
			virtual void GenerateData();
			virtual void FindMinMax();
			virtual void AllocateDisplacementField();

			//Protected Variables (In a Public Herence are Protected in the Son Class)
			int 															NumberOfParameters;
			int 															NumberOfBins;
			float															NoiseStdDev;
			bool															ban_send_vec;
			bool															ban_minmax;
			float															minimum;	
			float															maximum;
			float															minimumFixed;	
			float															maximumFixed;
			float															minimumMoving;	
			float															maximumMoving;

			IndexType													SamplesStepSize;
			SizeType														FixedImageSize;
			SizeType														MovingImageSize;
			VectorType 													ParticleParameters;	
			VectorType													ParticlePerturbations;
	
			typename VectorImageType::Pointer					DisplacementField;
			typename InternalImageType::SpacingType			MovingSpacing;
			typename InternalImageType::SpacingType			FixedSpacing;
			typename OutputImageType::Pointer					OutputImage;

		private:
			//Private Variables
			float															SimilarityValue;
			float															Weight;
			float 														c1,c2,c3;	//Auxiliar constants to Likelihhod
			bool 															m_PhysicalCoordinates;

			MatrixType 													AffineMatrix;

			typename InternalImageType::ConstPointer			MovingImage;
			typename InternalImageType::ConstPointer			FixedImage;
			typename InternalImageType::Pointer					AuxMovingImage;
			typename InputImageType::ConstPointer				InputFixedImage;
			typename InputImageType::ConstPointer				InputMovingImage;
			typename InterpolatorType::Pointer					m_Interpolator;

			//Private Functions
			void GenerateTransformationMatrix();
			void GenerateDisplacementField3D();
			void GenerateDisplacementField2D();
			float ComputeSimilarityMetric(int **Hist,int samples);
			void ComputeParticleWeight2D();
			void ComputeParticleWeight3D();

	 };//Class Particle


	//Resampling Node
	class ResamplingParticle{
		public:
			typedef itk::VariableLengthVector<float> VectorType;
			VectorType 			Parameters;		//Parameters Vector
			float 				Weight;
			float					Acc;				//Accumulate Density Function PDA
			int 					NumberOfParameters;

			//Builder
			ResamplingParticle(int number)
			{
				NumberOfParameters = number;
				Parameters.SetSize( NumberOfParameters );
			};

			~ResamplingParticle(){};//Destroy

			void SetParameters(VectorType Vec)
					{
						for(int i = 0; i < NumberOfParameters; i++)
							Parameters[i] = Vec[i];
					};
			//Return Parameters
			VectorType GetParameters()
							{
								return Parameters;
							};
	};

} //namespace ITK

#ifndef ITK_MANUAL_INSTANTIATION
#include "itkParticle.txx"
#endif

#endif

