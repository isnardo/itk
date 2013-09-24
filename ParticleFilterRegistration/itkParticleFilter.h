#ifndef _itkParticleFilter_h
#define _itkParticleFilter_h

#include <list>
#include "itkParticle.h"
#include "NumberOfCores.h"
//Multicore Processing
#include "itkMultiThreader.h"
#include "itkSimpleFastMutexLock.h"



namespace itk
{

	//Template Class Particle Filter
	template< class TInputImage, class TOutputImage >
	class ITK_EXPORT ParticleFilter : public Particle< TInputImage, TOutputImage >
	{
		public:
			// Standard class typedefs
			typedef ParticleFilter       												Self;
			typedef SmartPointer< Self >        									Pointer;
			typedef SmartPointer< const Self >        							ConstPointer;

			itkNewMacro(Self);
			itkTypeMacro(ParticleFilter, Particle);

			//Types definition
			typedef TInputImage 															InputImageType;
			typedef TOutputImage 														OutputImageType;

			//Const Macro to Define the Images Dimensions
			itkStaticConstMacro(ImageDimension, unsigned int,InputImageType::ImageDimension);

			typedef itk::Image< float,itkGetStaticConstMacro(ImageDimension)>		
																								InternalImageType;

 			typedef itk::CastImageFilter<InputImageType,InternalImageType >  																									InputCasterType;
 			typedef itk::CastImageFilter<InternalImageType,OutputImageType>  																									OutputCasterType;

			void SetFixedImage (const  InputImageType  *fixedImage  );
			void SetMovingImage(const  InputImageType *movingImage );

			virtual void SetInternalFixedImage (const  InternalImageType  *fixedImage )
										{
											FixedImage = fixedImage;
											this->FixedImageSize = FixedImage->GetLargestPossibleRegion().GetSize();
											this->FixedSpacing =  FixedImage->GetSpacing();
										}
			virtual void SetInternalMovingImage(const  InternalImageType *movingImage )
										{
											MovingImage = movingImage;
											this->MovingImageSize = MovingImage->GetLargestPossibleRegion().GetSize();
											this->MovingSpacing = MovingImage->GetSpacing();
										}

			typedef Particle<InternalImageType,InternalImageType>				ParticleType;

			typedef double 																CoordRepType;
			typedef InterpolateImageFunction<InternalImageType,CoordRepType> 	
																								InterpolatorType;
			typedef LinearInterpolateImageFunction<InternalImageType,CoordRepType >
																								DefaultInterpolatorType;
			//Particle Pointer
			typedef typename ParticleType::Pointer									ParticlePointer;
			typedef itk::Vector< float,itkGetStaticConstMacro(ImageDimension)>	
																								VectorPixelType;
			typedef itk::Image< VectorPixelType,
										itkGetStaticConstMacro(ImageDimension)>								
																								VectorImageType;

			//Threader Multicore			
			typedef itk::MultiThreader  												ThreaderType;


			itkSetObjectMacro( Interpolator, InterpolatorType );
			itkSetMacro( PhysicalCoordinates, bool);
			itkBooleanMacro( PhysicalCoordinates );

		
			//Set Number of Parameters to Estimated 2D(5 and 7) 3D(9 and 15)
			virtual void SetNumberOfParameters(int NumberOfParameters_t);
			//Set Number of Particles
			virtual void SetNumberOfParticles(int NumberOfParticles_t)
								{
									NumberOfParticles = NumberOfParticles_t;
								};
			//Set Maximum Number of Iterations
			virtual void SetMaxIterations(int MaxIterations_t)
								{
									MaxIterations = MaxIterations_t;
								};
			virtual void SetThresholdStop(float ThresholdStop_t)
								{
									ThresholdStop = ThresholdStop_t;
								};
			//Return Transformed Image
			typename OutputImageType::Pointer GetTransformedImage();
			//Return Estimated Parameters
			typename ParticleType::VectorType GetParameters()
								{
									return Estimated->GetParameters();
								};

			virtual typename VectorImageType::Pointer GetDisplacementField(){
																		Estimated->GenerateDisplacementField();
																		return Estimated->GetDisplacementField();
																	};
			virtual void SetIntialParameters(typename ParticleType::VectorType Parameters_t)
									{
										this->SetParameters( Parameters_t );
										flag_parameters = true;
									};

		protected:
			ParticleFilter();
			~ParticleFilter()
				{
					ParticleList.clear();
					ResamplingList.clear();
				};

			//Main Function of the Class
			virtual void GenerateData();

			//Protected Variables
			int 																	NumberOfParticles;
			int 																	MaxIterations;
			int 																	ThresholdStop;	
			bool																	flag_parameters;

			ParticlePointer													Estimated;

		private:
			//Private Variables
			float 																IntensitiesRange;
			bool 																	m_PhysicalCoordinates;

			void ComputeIntensitiesRange();
			void AllocateParticlesList();
			float ComputePerturbationsDegeneration();

			typename InterpolatorType::Pointer							m_Interpolator;
			typename InputImageType::ConstPointer						InputFixedImage;
			typename InputImageType::ConstPointer						InputMovingImage;

			//List to construct a lot of particles
			std::list<ParticlePointer>										ParticleList;
			typename std::list<ParticlePointer>::iterator			ParticleIterator;

			//List to Resample Particles
			std::list<ResamplingParticle*>								ResamplingList;
			typename std::list<ResamplingParticle*>::iterator		ResamplingIterator;


			typename InternalImageType::ConstPointer					MovingImage;
			typename InternalImageType::ConstPointer					FixedImage;

			void ComputeWeigthsMultiCore();

			static ITK_THREAD_RETURN_TYPE ThreaderCallback( void * arg );
	};//Class Particle Filter


/*
	static itk::SimpleFastMutexLock mutex;
	//Class to implement MultiCore Particle Filtering
	class MultiCoreParticleStructure{

		public:
			void*			Iterator;
			void*			IteratorEnd;

			MultiCoreParticleStructure(){};
			~MultiCoreParticleStructure(){};

	};

//MULTICORE
	const static int GDimension = 3;
	typedef itk::Image< float,GDimension>												GImageType;
	typedef Particle<GImageType,GImageType>											GParticleType;
	typedef typename GParticleType::Pointer											GParticlePointer;
	std::list<GParticlePointer>															GParticleList;
	typename std::list<GParticlePointer>::iterator									GParticleIterator;


	int Gcont = 0;
	static itk::SimpleFastMutexLock mutex;
	ITK_THREAD_RETURN_TYPE ThreaderCallback( void * arg );
*/
	
} //namespace ITK

#ifndef ITK_MANUAL_INSTANTIATION
#include "itkParticleFilter.txx"
#endif

#endif
