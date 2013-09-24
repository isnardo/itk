#ifndef _itkParticleFilter_txx
#define _itkParticleFilter_txx

#include "itkParticleFilter.h"

namespace itk
{
	//**************************************
	//Functions definition of Particle Class
	//**************************************
		
	//Class Builder
	template< class TInputImage, class TOutputImage>	
	ParticleFilter< TInputImage, TOutputImage>
	::ParticleFilter()
	{
		this->SetNumberOfRequiredInputs(0);//Number of requaired images to Update

		//Noise Standar Deviation
		this->NoiseStdDev = 0.2;		
		this->SetNoiseStandardDeviation(this->NoiseStdDev);

		this->ban_send_vec = false;			//Band to Send Vector Field
		this->NumberOfBins = 32;				//Number of Bins

		NumberOfParticles = 500;
		MaxIterations 		= 600;
		ThresholdStop 		= 1e-5;	

		//Default Internal Interpolator
		typename DefaultInterpolatorType::Pointer   Interp = DefaultInterpolatorType::New();
	  	m_Interpolator =  static_cast<InterpolatorType*>( Interp.GetPointer() );

		//StepSize and Parameters
		if( this->ImageDimension == 3)
		{
			this->NumberOfParameters = 9;	//Number of Parameteres
			this->SamplesStepSize[0] = 8;
			this->SamplesStepSize[1] = 8;
			this->SamplesStepSize[2] = 3;
		}
		else
		{
			this->NumberOfParameters = 5;	//Number of Parameteres
			this->SamplesStepSize[0] = 8;
			this->SamplesStepSize[1] = 8;
		}

		//Inicialization of Perturbations
		SetNumberOfParameters( this->NumberOfParameters );
		flag_parameters = false;
		m_PhysicalCoordinates = false;
	}


	//Set Fixed Image (Input 0)
	template< class TInputImage, class TOutputImage>	
	void ParticleFilter< TInputImage, TOutputImage>
	::SetFixedImage(const InputImageType* fixedImage)
	{
		//Fixed Image (Input 0)
		//this->SetNthInput(0, const_cast<TInputImage*>(fixedImage));

		InputFixedImage = fixedImage;

		//Cast Fixed Imgae To Internal ImageType Float
		typename InputCasterType::Pointer ImageCaster = InputCasterType::New();
		ImageCaster->SetInput( fixedImage );
		//Update Caster
		try
		{
			ImageCaster->Update();
		}
		catch( itk::ExceptionObject & excep )
		{
			std::cerr << "Exception catched !" << std::endl;
			std::cerr << excep << std::endl;
		}

		FixedImage = ImageCaster->GetOutput();
		this->FixedImageSize = FixedImage->GetLargestPossibleRegion().GetSize();
		this->FixedSpacing = FixedImage->GetSpacing();
	}

	//Set Moving Image (Input 1)
	template< class TInputImage, class TOutputImage>	
	void ParticleFilter< TInputImage, TOutputImage>
	::SetMovingImage(const TInputImage* movingImage)
	{
		//Moving Image (Input 1)
		//this->SetNthInput(1, const_cast<TInputImage*>(movingImage));
		InputMovingImage = movingImage;

		//Cast Moving Imgae To Internal ImageType Float
		typename InputCasterType::Pointer ImageCaster = InputCasterType::New();
		ImageCaster->SetInput( movingImage );
		//Update Caster
		try
		{
			ImageCaster->Update();
		}
		catch( itk::ExceptionObject & excep )
		{
			std::cerr << "Exception catched !" << std::endl;
			std::cerr << excep << std::endl;
		}
		MovingImage = ImageCaster->GetOutput();
		this->MovingImageSize = MovingImage->GetLargestPossibleRegion().GetSize();
		this->MovingSpacing = MovingImage->GetSpacing();
	}

	//Compute Intesities Range
	template< class TInputImage, class TOutputImage>	
	void ParticleFilter< TInputImage, TOutputImage>
	::ComputeIntensitiesRange()
	{
		//Compute Intesities Range With inherited Function
		this->FindMinMax();
		IntensitiesRange = this->maximum - this->minimum;
	}

	//Set Number of Parameter and allocate ParticlePerturbations
	template< class TInputImage, class TOutputImage>	
	void ParticleFilter< TInputImage, TOutputImage>
	::SetNumberOfParameters(int NumberOfParameters_t)
	{
		// 5 Parameters  - Affine Without shearing   2D
 		// 7 Parameters  - Affine Wit shearing       2D
		// 9 Parameters  - Affine Without shearing	3D
		// 15 Parameters - Affine With shearing	   3D

		this->NumberOfParameters = NumberOfParameters_t;
		this->ParticlePerturbations.SetSize( this->NumberOfParameters );

		//Volume
		if ( this->ImageDimension == 3 )
		{
			this->ParticlePerturbations[0] = 1.0;	//Angle x in Degrees
			this->ParticlePerturbations[1] = 1.0;	//Angle y in Degrees
			this->ParticlePerturbations[2] = 1.0;	//Angle z in Degrees

			this->ParticlePerturbations[3] = 10.0;	//X Shift
			this->ParticlePerturbations[4] = 10.0;	//Y Shift
			this->ParticlePerturbations[5] = 1.0;	//Z Shift

			this->ParticlePerturbations[6] = 2.5e-3;	//X Scale
			this->ParticlePerturbations[7] = 2.5e-3;	//Y Scale  
			this->ParticlePerturbations[8] = 2.5e-3;	//Z Scale

			if( this->NumberOfParameters == 15 ){ //With Shearing
				this->ParticlePerturbations[9] = 2.5e-5; 	//XY Shearing
				this->ParticlePerturbations[10] = 2.5e-5; //XZ Shearing
				this->ParticlePerturbations[11] = 2.5e-5; //YX Shearing
				this->ParticlePerturbations[12] = 2.5e-5; //YZ Shearing
				this->ParticlePerturbations[13] = 2.5e-5; //ZX Shearing
				this->ParticlePerturbations[14] = 2.5e-5; //ZY Shearing
			}
		}
		else
		{					
			this->ParticlePerturbations[0] = 15.0;	//Angle in Degrees
			this->ParticlePerturbations[1] = 10.0;	//X Shift
			this->ParticlePerturbations[2] = 10.0;	//Y Shift
			this->ParticlePerturbations[3] = 0.0025;	//X Scale
			this->ParticlePerturbations[4] = 0.0025;	//Y Scale

			if( this->NumberOfParameters == 7)//With Shearing
			{
				this->ParticlePerturbations[5] = 2.5e-5;	//X Shearing
				this->ParticlePerturbations[6] = 2.5e-5;	//Y Shearing
			}
		}
	}

	//Allocate Particles List
	template< class TInputImage, class TOutputImage>	
	void ParticleFilter< TInputImage, TOutputImage>
	::AllocateParticlesList()
	{
		//Pointers
		ParticlePointer  			p;
		ResamplingParticle		*r;
		

		//Generate Particles
		for(int i = 0; i < NumberOfParticles; i++){

			r = new ResamplingParticle( this->NumberOfParameters );

			p = ParticleType::New();	//New particle
			//Set Values to each particle
			p->SetInternalFixedImage( FixedImage );
			p->SetInternalMovingImage( MovingImage );
			p->SetSamplesStepSize( this->SamplesStepSize );
			p->SetNumberOfBins( this->NumberOfBins );
			p->SetNoiseStandardDeviation( this->NoiseStdDev );
			p->SetNumberOfParameters( this->NumberOfParameters );
			if( m_PhysicalCoordinates )
				p->PhysicalCoordinatesOn();

		//	p->SetInterpolator( static_cast<InterpolatorType*>
		//							 (m_Interpolator->CreateAnother().GetPointer()) );
			p->SetInterpolator( m_Interpolator );

			if( flag_parameters )
				p->SetParameters( this->ParticleParameters );

			//Add each element to the Lists
			ParticleList.push_back( p );
			ResamplingList.push_back( r );
		}

		//Set values to the estimated
		Estimated = ParticleType::New();
		Estimated->SetInternalFixedImage( FixedImage );
		Estimated->SetInternalMovingImage( MovingImage );
		Estimated->SetNumberOfParameters( this->NumberOfParameters );
		Estimated->SetSamplesStepSize( this->SamplesStepSize );
		Estimated->SetNumberOfBins( this->NumberOfBins );
		Estimated->SetNoiseStandardDeviation( this->NoiseStdDev );
		Estimated->SetInterpolator( m_Interpolator );

	}


	//Compute Perturbations Degeneration
	template< class TInputImage, class TOutputImage>	
	float ParticleFilter< TInputImage, TOutputImage>
	::ComputePerturbationsDegeneration()
	{
		float d = 0.0;
		for(int i = 0; i < this->NumberOfParameters; i++)
			d += this->ParticlePerturbations[i]*this->ParticlePerturbations[i]; 
		d = sqrt(d);
		return d;
	}

	//Return the Transformed Image
	template< class TInputImage, class TOutputImage>	
	typename TOutputImage::Pointer ParticleFilter< TInputImage, TOutputImage>
	::GetTransformedImage()
	{
		//Cast To Output ImgeT ype
		typename OutputCasterType::Pointer ImageCaster = OutputCasterType::New();
		ImageCaster->SetInput( Estimated->GetTransformedImage() );
		//Update Caster
		try
		{
			ImageCaster->Update();
		}
		catch( itk::ExceptionObject & excep )
		{
			std::cerr << "Exception catched !" << std::endl;
			std::cerr << excep << std::endl;
		}

		return ImageCaster->GetOutput();
	}

	//********************************************
	//Generate Data   (Main fucntion of the Class)
	//********************************************
	template< class TInputImage, class TOutputImage>	
	void ParticleFilter< TInputImage, TOutputImage>
	::GenerateData()
	{
		float sum;
		float w;
		float Neff,aux;
		float facDeg = exp(log(1e-5)/MaxIterations);
		float randomNumberU;

		srand ( time(NULL) );			//Start Random Seed

		//Allocate Particle List
		AllocateParticlesList();

		std::cout<<" Number of Particles: "<< this->NumberOfParticles <<std::endl;  
		std::cout<<" Number of Iterations: "<< this->MaxIterations <<std::endl;  
		std::cout<<" Number of Parameters: "<< this->NumberOfParameters<<std::endl;
		std::cout<<" Physical Coordinates:";
		if( m_PhysicalCoordinates )
			std::cout<<" On"<<std::endl;
		else
			std::cout<<" Off (Computing Voxel coordinates)"<<std::endl;
		std::cout<<" Histograms' bins: "<< this->NumberOfBins<<std::endl;
		std::cout<<" Sampling Spacing: "<< this->SamplesStepSize <<std::endl;
		std::cout<<" Initial perturbations' variances: "<< this->ParticlePerturbations <<std::endl;
		
		//Start Particle Filtering
		for(int ite = 0; ite < this->MaxIterations; ite++)
		{

			sum = 0.0;
			//Perturb Particles
			for (ParticleIterator = ParticleList.begin(); 
					ParticleIterator != ParticleList.end(); ParticleIterator++)
			{
				//Perturb Particles		
				(*ParticleIterator)->UpsetParticle( this->ParticlePerturbations );	
				//Compute Weigth
				w = (*ParticleIterator)->ComputeParticleWeight();	
				sum += w; //Add weights
			}

/*
//MULTICORE
			ComputeWeigthsMultiCore();

			sum = 0.0;
			//Compute Sum of Weights
			for (ParticleIterator = ParticleList.begin(); 
					ParticleIterator != ParticleList.end(); ParticleIterator++)
			{
				//Compute Particle Weigth
				w = (*ParticleIterator)->GetParticleWeight();	
				sum += w; //Add weights
			}

*/
	
			//Weights Sum for Normalization
			sum = 1.0/sum;

			//Computing PDA
			for (ParticleIterator = ParticleList.begin(), ResamplingIterator = ResamplingList.begin(); 
					ParticleIterator != ParticleList.end(); 
					ParticleIterator++,ResamplingIterator++)
			{
				//Normalizing weights
				(*ResamplingIterator)->Weight = (*ParticleIterator)->GetParticleWeight() * sum;
				//Making a copy of parameters for resampling
				(*ResamplingIterator)->SetParameters( (*ParticleIterator)->GetParameters() );

				//Obtain Accumulate Density Function
				if(ParticleIterator == ParticleList.begin())
					(*ResamplingIterator)->Acc = (*ResamplingIterator)->Weight;
				else
					(*ResamplingIterator)->Acc = w + (*ResamplingIterator)->Weight;

				//Previous PDA
				w = (*ResamplingIterator)->Acc;		
			}

			Neff = 0.0;

			//Resampling Particles
			for(ParticleIterator = ParticleList.begin(); 
					ParticleIterator != ParticleList.end(); 
					ParticleIterator++)
			{
				//Random number with Uniform distribution [0,1]
				randomNumberU = rand_uniform();
				//Search into PDA the new particle	
				for ( ResamplingIterator = ResamplingList.begin(); 
						ResamplingIterator != ResamplingList.end(); ResamplingIterator++)
				{
					if( (*ResamplingIterator)->Acc >= randomNumberU )
					{
						(*ParticleIterator)->SetParameters( (*ResamplingIterator)->Parameters );
						(*ParticleIterator)->SetParticleWeight( (*ResamplingIterator)->Weight );
						break;
					}
				}
			}//End Particles Resampling

			Neff = ComputePerturbationsDegeneration();

			if( ite == 0 )
				this->ReducePerturbations( 0.1 );
			else
				this->ReducePerturbations( facDeg );

			Estimated->FillParameters(0.0);
			for(ParticleIterator = ParticleList.begin(); 
				ParticleIterator != ParticleList.end(); ParticleIterator++)
			{
				Estimated->ScaleAddVector( (*ParticleIterator)->GetParameters(),
													(*ParticleIterator)->GetParticleWeight() );
			}

			//Degeneration perturbations criteria
			if( Neff < ThresholdStop )
				break;

			if( ( ite+1 )%100 == 0)
			{
				std::cout<<" Iteration: "<<ite+1<<std::endl;
				std::cout<<" Deg. Factor: "<<Neff<<std::endl;
			   std::cout<<" Estimated: "<<GetParameters()<<std::endl;
			}
		}


	}	

/*
//MULTICORE
	//Weigths Process Multicore
	template< class TInputImage, class TOutputImage>	
	void ParticleFilter< TInputImage, TOutputImage>
	::ComputeWeigthsMultiCore()
	{
		ThreaderType::Pointer threader = ThreaderType::New();
		const unsigned int numberOfCores = getNumberOfCores();
	   const unsigned int numberOfThreads = numberOfCores;
		typename std::list<ParticlePointer>::iterator			IteratorEnd;

		MultiCoreParticleStructure Datos;

		ParticleIterator = ParticleList.begin();	
		ParticleIterator--;	
		IteratorEnd = ParticleList.end();	

		Datos.Iterator = (void*)&ParticleIterator;
		Datos.IteratorEnd = (void*)&IteratorEnd;

		threader->SetNumberOfThreads( numberOfThreads );
   	threader->SetSingleMethod( this->ThreaderCallback, &Datos);

   	threader->SingleMethodExecute();
	}


//MULTICORE
	template< class TInputImage, class TOutputImage>	
	ITK_THREAD_RETURN_TYPE ParticleFilter< TInputImage, TOutputImage>
	::ThreaderCallback( void * arg )
	{

		typedef typename std::list<ParticlePointer>::iterator		IteratorType;
		IteratorType*			Iterator;
		IteratorType*			IteratorEnd;
		MultiCoreParticleStructure* Datos = static_cast< MultiCoreParticleStructure * >( arg );
	
		Iterator = static_cast<  IteratorType* >(Datos->Iterator);
		IteratorEnd = static_cast<  IteratorType* >(Datos->IteratorEnd);


		while( 1 )
		{
//			mutex.Lock();
			++Iterator;
//			mutex.Unlock();
			std::cout<<Iterator<<std::endl;
/*			if( Iterator != IteratorEnd )
				(*(*Iterator))->ComputeParticleWeight();	
			else 
				break;

		}


		ITK_THREAD_RETURN_TYPE value;
		return value;
	}*/

}//namespace itk 

#endif
