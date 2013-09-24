#ifndef _itkParticle_txx
#define _itkParticle_txx

#include "itkParticle.h"

namespace itk
{
	//**************************************
	//Functions definition of Particle Class
	//**************************************
		
	//Class Builder
	template< class TInputImage, class TOutputImage>	
	Particle< TInputImage, TOutputImage>
	::Particle()
	{
		this->SetNumberOfRequiredInputs(0);//Number of requaired images to Update

		////****** Default Values *******///

		//Noise Standar Deviation
		NoiseStdDev = 0.2;		
		SetNoiseStandardDeviation(NoiseStdDev);

		ban_send_vec = false;			//Flag to Send Vector Field
		ban_minmax   = false;			//Flag to min max computation
		NumberOfBins = 256;				//Number of Bins
		Weight = 0.0;
		m_PhysicalCoordinates = false;

		//Default Interpolator
		typename DefaultInterpolatorType::Pointer   Interp = DefaultInterpolatorType::New();
	  	m_Interpolator =  static_cast<InterpolatorType*>( Interp.GetPointer() );

		//StepSize and Parameters
		if( ImageDimension == 3)
		{
			NumberOfParameters = 9;	//Number of Parameteres
			SamplesStepSize[0] = 8;
			SamplesStepSize[1] = 8;
			SamplesStepSize[2] = 4;
		}
		else
		{
			NumberOfParameters = 5;	//Number of Parameteres
			SamplesStepSize[0] = 8;
			SamplesStepSize[1] = 8;
		}


	}

	//Set Fixed Image (Input 0)
	template< class TInputImage, class TOutputImage>	
	void Particle< TInputImage, TOutputImage>
	::SetFixedImage(const InputImageType* fixedImage)
	{
		InputFixedImage = fixedImage;

		//Cast Fixed Imgae To InternalImageType Float
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
		FixedImageSize = FixedImage->GetLargestPossibleRegion().GetSize();
		FixedSpacing =  FixedImage->GetSpacing();
	}

	//Set Moving Image (Input 1)
	template< class TInputImage, class TOutputImage>	
	void Particle< TInputImage, TOutputImage>
	::SetMovingImage(const TInputImage* movingImage)
	{
		InputMovingImage = movingImage;

		//Cast Moving Imgae To InternalImageType Float
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
		MovingImageSize = MovingImage->GetLargestPossibleRegion().GetSize();
		MovingSpacing =  MovingImage->GetSpacing();
	}


	//Set Moving Image (Input 1)
	template< class TInputImage, class TOutputImage>	
	void Particle< TInputImage, TOutputImage>
	::SetNoiseStandardDeviation(float StdDev)
	{
		NoiseStdDev = StdDev;
		c1 = 1.0/(sqrt(2.0*M_PI)*NoiseStdDev);
		c2 = -1.0/(2.0*NoiseStdDev*NoiseStdDev);
	}

	//Set Parameters Vector
	template< class TInputImage, class TOutputImage>	
	void Particle< TInputImage, TOutputImage>
	::SetParameters(VectorType Vector)
	{
		ParticleParameters.SetSize( NumberOfParameters );

		for(int i = 0; i < NumberOfParameters; i++)
			ParticleParameters[i] = Vector[i];
	}

	//Set Variance of Perturbations
	template< class TInputImage, class TOutputImage>	
	void Particle< TInputImage, TOutputImage>
	::SetVarianceOfPerturbations(VectorType Vector)
	{

		ParticlePerturbations.SetSize( NumberOfParameters );

		for(int i = 0; i < NumberOfParameters; i++)
			ParticlePerturbations[i] = Vector[i];
	}

	//Set Parameters Vector
	template< class TInputImage, class TOutputImage>	
	void Particle< TInputImage, TOutputImage>
	::SetNumberOfParameters(int NumberOfParameters_t)
	{
		// 5 Parameters  - Affine Without shearing   2D
 		// 7 Parameters  - Affine Wit shearing       2D
		// 9 Parameters  - Affine Without shearing	3D
		// 15 Parameters - Affine With shearing	   3D

		NumberOfParameters = NumberOfParameters_t;

		ParticleParameters.SetSize( NumberOfParameters );

		//Volume
		if ( ImageDimension == 3)
		{
			ParticleParameters[0] = 0.0;	//Angle x in Degrees
			ParticleParameters[1] = 0.0;	//Angle y in Degrees
			ParticleParameters[2] = 0.0;	//Angle z in Degrees

			ParticleParameters[3] = 0.0; 	//X Shift
			ParticleParameters[4] = 0.0;	//Y Shift
			ParticleParameters[5] = 0.0;	//Z Shift

			ParticleParameters[6] = 1.0;	//X Scale
			ParticleParameters[7] = 1.0;	//Y Scale  
			ParticleParameters[8] = 1.0;	//Z Scale

			if( NumberOfParameters == 15 ){ //With Shearing
				ParticleParameters[9] = 0.0; 	//XY Shearing
				ParticleParameters[10] = 0.0; //XZ Shearing
				ParticleParameters[11] = 0.0; //YX Shearing
				ParticleParameters[12] = 0.0; //YZ Shearing
				ParticleParameters[13] = 0.0; //ZX Shearing
				ParticleParameters[14] = 0.0; //ZY Shearing
			}
		}
		//2D
		else
		{					
			ParticleParameters[0] = 0.0;	//Angle in Degrees
			ParticleParameters[1] = 0.0;	//X Shift
			ParticleParameters[2] = 0.0;	//Y Shift
			ParticleParameters[3] = 1.0; 	//X Scale
			ParticleParameters[4] = 1.0;	//Y Scale

			if( NumberOfParameters == 7)//With Shearing
			{
				ParticleParameters[5] = 0.0;	//X Shearing
				ParticleParameters[6] = 0.0;	//Y Shearing
			}
		}
	}

	//Generate Transformation Matrix
	template< class TInputImage, class TOutputImage>	
	void Particle< TInputImage, TOutputImage>
	::GenerateTransformationMatrix()
	{
		float cosX,sinX,cosY,sinY,cosZ,sinZ,cosT,sinT;		//sinc functions
		float thetaX_rad, thetaY_rad, thetaZ_rad;	//angles in radians
		float dx,dy,dz;						//Displacements
		float ex,ey,ez;						//Scales
		float hxy,hxz,hyx,hyz,hzx,hzy,hx,hy;	//shearing
		float ang = M_PI/180.0;

		//Compute the Transformation Matrix
		switch( NumberOfParameters )
		{
			case 9: //Volume Affine without shearing

				//Angles in radians
				thetaX_rad = ParticleParameters[0]*ang;  //X angle
				thetaY_rad = ParticleParameters[1]*ang;  //Y angle
				thetaZ_rad = ParticleParameters[2]*ang;  //Z angle
				//cos and sin
				cosX = cos(thetaX_rad);		
				sinX = sin(thetaX_rad);		
				cosY = cos(thetaY_rad);		
				sinY = sin(thetaY_rad);		
				cosZ = cos(thetaZ_rad);		
				sinZ = sin(thetaZ_rad);		
				//Displacements
				dx = ParticleParameters[3];
				dy = ParticleParameters[4];
				dz = ParticleParameters[5];
				//Scales
				ex = ParticleParameters[6];
				ey = ParticleParameters[7];
				ez = ParticleParameters[8];

				//Generate transformation Matrix
				AffineMatrix(0,0) = ex*cosY*cosZ;
				AffineMatrix(0,1) = -ey*(cosX*sinZ - cosZ*sinX*sinY);
				AffineMatrix(0,2) = ez*(sinX*sinZ + cosX*cosZ*sinY);
				AffineMatrix(0,3) = dx;

				AffineMatrix(1,0) = ex*cosY*sinZ;
				AffineMatrix(1,1) = ey*(cosX*cosZ + sinX*sinY*sinZ);
				AffineMatrix(1,2) = -ez*(cosZ*sinX - cosX*sinY*sinZ);
				AffineMatrix(1,3) = dy;

				AffineMatrix(2,0) = -ex*sinY;
				AffineMatrix(2,1) = ey*cosY*sinX;
				AffineMatrix(2,2) = ez*cosX*cosY;
				AffineMatrix(2,3) = dz;

				AffineMatrix(3,0) = 0.0;
				AffineMatrix(3,1) = 0.0;
				AffineMatrix(3,2) = 0.0;
				AffineMatrix(3,3) = 1.0;

				break;

			case 5: //2D Affine without shearing

				//Angle in radians
				thetaX_rad = ParticleParameters[0]*ang;
				//cos and sin
				cosT = cos(thetaX_rad);		
				sinT = sin(thetaX_rad);		
				//Displacements
				dx = ParticleParameters[1];
				dy = ParticleParameters[2];
				//Scales
				ex = ParticleParameters[3];
				ey = ParticleParameters[4];

				//Generate transformation Matrix
				AffineMatrix(0,0) = ex*cosT;
				AffineMatrix(0,1) = -ey*sinT;
				AffineMatrix(0,2) = dx;

				AffineMatrix(1,0) = ex*sinT;
				AffineMatrix(1,1) = ey*cosT;
				AffineMatrix(1,2) = dy;

				AffineMatrix(2,0) = 0.0;
				AffineMatrix(2,1) = 0.0;
				AffineMatrix(2,2) = 1.0;
			
				break;

			case 15: //Volume Affine with shearing

				//Angles in radians
				thetaX_rad = ParticleParameters[0]*ang;  //X angle
				thetaY_rad = ParticleParameters[1]*ang;  //Y angle
				thetaZ_rad = ParticleParameters[2]*ang;  //Z angle
				//cos and sin
				cosX = cos(thetaX_rad);		
				sinX = sin(thetaX_rad);		
				cosY = cos(thetaY_rad);		
				sinY = sin(thetaY_rad);		
				cosZ = cos(thetaZ_rad);		
				sinZ = sin(thetaZ_rad);		
				//Displacements
				dx = ParticleParameters[3];
				dy = ParticleParameters[4];
				dz = ParticleParameters[5];
				//Scales
				ex = ParticleParameters[6];
				ey = ParticleParameters[7];
				ez = ParticleParameters[8];

				//Shearing
				hxy = ParticleParameters[9];
				hxz = ParticleParameters[10];
				hyx = ParticleParameters[11];
				hyz = ParticleParameters[12];
				hzx = ParticleParameters[13];
				hzy = ParticleParameters[14];

				//Generate transformation Matrix
				AffineMatrix(0,0) = ex*(cosY*cosZ - hyx*(cosX*sinZ - cosZ*sinX*sinY) 
																+ hzx*(sinX*sinZ + cosX*cosZ*sinY));
				AffineMatrix(0,1) = ey*(hzy*(sinX*sinZ + cosX*cosZ*sinY) - cosX*sinZ 
																+ hxy*cosY*cosZ + cosZ*sinX*sinY);
				AffineMatrix(0,2) = ez*(sinX*sinZ - hyz*(cosX*sinZ - cosZ*sinX*sinY) 
																+ hxz*cosY*cosZ + cosX*cosZ*sinY);
				AffineMatrix(0,3) = dx;

				AffineMatrix(1,0) = ex*(cosY*sinZ + hyx*(cosX*cosZ + sinX*sinY*sinZ) 
																- hzx*(cosZ*sinX - cosX*sinY*sinZ));
				AffineMatrix(1,1) = ey*(cosX*cosZ - hzy*(cosZ*sinX - cosX*sinY*sinZ) 
																+ sinX*sinY*sinZ + hxy*cosY*sinZ);
				AffineMatrix(1,2) = ez*(hyz*(cosX*cosZ + sinX*sinY*sinZ) - cosZ*sinX 
																+ hxz*cosY*sinZ + cosX*sinY*sinZ);
				AffineMatrix(1,3) = dy;

				AffineMatrix(2,0) = ex*(hzx*cosX*cosY - sinY + hyx*cosY*sinX);
				AffineMatrix(2,1) = ey*(cosY*sinX - hxy*sinY + hzy*cosX*cosY);
				AffineMatrix(2,2) = ez*(cosX*cosY - hxz*sinY + hyz*cosY*sinX);
				AffineMatrix(2,3) = dz;

				AffineMatrix(3,0) = 0.0;
				AffineMatrix(3,1) = 0.0;
				AffineMatrix(3,2) = 0.0;
				AffineMatrix(3,3) = 1.0;

				break;

			case 7: //2D Affine with shearing

				//Angle in radians
				thetaX_rad = ParticleParameters[0]*ang;
				//cos and sin
				cosT = cos(thetaX_rad);		
				sinT = sin(thetaX_rad);		
				//Displacements
				dx = ParticleParameters[1];
				dy = ParticleParameters[2];
				//Scales
				ex = ParticleParameters[3];
				ey = ParticleParameters[4];
				//Shearing
				hx = ParticleParameters[5];
				hy = ParticleParameters[6];

				//Generate transformation Matrix
				AffineMatrix(0,0) = ex*(cosT - hy*sinT);
				AffineMatrix(0,1) = -ey*(sinT - hx*cosT);
				AffineMatrix(0,2) = dx;

				AffineMatrix(1,0) = ex*(sinT + hy*cosT);
				AffineMatrix(1,1) = ey*(cosT + hx*sinT);
				AffineMatrix(1,2) = dy;

				AffineMatrix(2,0) = 0.0;
				AffineMatrix(2,1) = 0.0;
				AffineMatrix(2,2) = 1.0;

				break;
			
		}//End Switch
	}

	//Reduce Perturbations
	template< class TInputImage, class TOutputImage>	
	void Particle< TInputImage, TOutputImage>
	::ReducePerturbations(float facDeg)
	{
		//Reduce variance of perturbations
		for(int i = 0; i < NumberOfParameters; i++)
			ParticlePerturbations[i] = ParticlePerturbations[i] * facDeg; 
	}

	//Upset randomly parameters vector of the particle 
	template< class TInputImage, class TOutputImage>	
	void Particle< TInputImage, TOutputImage>
	::UpsetParticle()
	{
		for(int i = 0; i < NumberOfParameters; i++)
			ParticleParameters[i] += rand_normal( ParticlePerturbations[i] ); 
	}
	
	//Upset randomly parameters vector of the particle with a given Perturbations Parameters
	template< class TInputImage, class TOutputImage>	
	void Particle< TInputImage, TOutputImage>
	::UpsetParticle(VectorType Perturbations)
	{
		for(int i = 0; i < NumberOfParameters; i++)
			ParticleParameters[i] += rand_normal( Perturbations[i] ); 
	}

	//Scale Vector and Add too the actual
	template< class TInputImage, class TOutputImage>	
	void Particle< TInputImage, TOutputImage>
	::ScaleAddVector(VectorType Vec, float fac)
	{
		//Scale and Add Vector
		for(int i = 0; i < NumberOfParameters; i++)
			ParticleParameters[i] += Vec[i] * fac; 
	}

	//Allocate Vector Field
	template< class TInputImage, class TOutputImage>	
	void Particle< TInputImage, TOutputImage>
	::AllocateDisplacementField()
	{
		//Allocate Memory for Displacementes Vetor Field
		DisplacementField = VectorImageType::New();
		DisplacementField->SetSpacing( FixedImage->GetSpacing() );
		DisplacementField->SetOrigin( FixedImage->GetOrigin() );
		DisplacementField->SetDirection( FixedImage->GetDirection() );
		DisplacementField->SetRegions( FixedImage->GetLargestPossibleRegion() );
		DisplacementField->Allocate();
	}

	//Generate Vector Field
	template< class TInputImage, class TOutputImage>	
	void Particle< TInputImage, TOutputImage>
	::GenerateDisplacementField(){
		//If Vector Field Memory was not reserved
		if( !DisplacementField )
			AllocateDisplacementField();

		//Compute Transformation Matrix
		GenerateTransformationMatrix();

		if( ImageDimension == 3)
			GenerateDisplacementField3D();
		else
			GenerateDisplacementField2D();

	}

	//Generate Vector Field 2D
	template< class TInputImage, class TOutputImage>	
	void Particle< TInputImage, TOutputImage>
	::GenerateDisplacementField2D(){

		float cx,cy;
		float A00,A01,A02;
		float A10,A11,A12;
		float ax,ay;

		VectorIndexType     Index;
		VectorPixelType     Shift;

		cx = MovingImageSize[0]/2.0;	//X-axis center
		cy = MovingImageSize[1]/2.0;	//Y-axis center

		A02 = AffineMatrix(0,2);	//Dx
		A12 = AffineMatrix(1,2);	//Dy
		
		for(Index[1] = 0; Index[1] < FixedImageSize[1]; Index[1]++)
		{
			ay = Index[1] - cy;
			A01 = AffineMatrix(0,1) * ay;
			A11 = AffineMatrix(1,1) * ay;

			for(Index[0] = 0; Index[0] < FixedImageSize[0]; Index[0]++)
			{
				ax = Index[0] - cx;
				A00 = AffineMatrix(0,0) * ax;
				A10 = AffineMatrix(1,0) * ax;

				Shift[0] = A00 + A01 + A02 + cx - Index[0];
				Shift[1] = A10 + A11 + A12 + cy - Index[1];
				
				DisplacementField->SetPixel( Index,Shift );
			}//End X-Axis
		}//End Y-Axis
	}

	//Generate Vector Field 3D
	template< class TInputImage, class TOutputImage>	
	void Particle< TInputImage, TOutputImage>
	::GenerateDisplacementField3D(){
		float cx,cy,cz;
		float A00,A01,A02,A03;
		float A10,A11,A12,A13;
		float A20,A21,A22,A23;
		float ax,ay,az;

		VectorIndexType     Index;
		VectorPixelType     Shift;

		cx = MovingImageSize[0]/2.0;	//X-axis center
		cy = MovingImageSize[1]/2.0;	//Y-axis center
		cz = MovingImageSize[2]/2.0;	//Z-axis center

		A03 = AffineMatrix(0,3);	//Dx
		A13 = AffineMatrix(1,3);	//Dy
		A23 = AffineMatrix(2,3);	//Dz

		for(Index[2] = 0; Index[2] < FixedImageSize[2]; Index[2]++)
		{
			az = Index[2] - cz;
			A02 = AffineMatrix(0,2) * az;
			A12 = AffineMatrix(1,2) * az;
			A22 = AffineMatrix(2,2) * az;

			for(Index[1] = 0; Index[1] < FixedImageSize[1]; Index[1]++)
			{
				ay = Index[1] - cy;
				A01 = AffineMatrix(0,1) * ay;
				A11 = AffineMatrix(1,1) * ay;
				A21 = AffineMatrix(2,1) * ay;

				for(Index[0] = 0; Index[0] < FixedImageSize[0]; Index[0]++)
				{
					ax = Index[0] - cx;
					A00 = AffineMatrix(0,0) * ax;
					A10 = AffineMatrix(1,0) * ax;
					A20 = AffineMatrix(2,0) * ax;

					Shift[0] = A00 + A01 + A02 + A03 + cx - Index[0];
					Shift[1] = A10 + A11 + A12 + A13 + cy - Index[1];
					Shift[2] = A20 + A21 + A22 + A23 + cz - Index[2];

					DisplacementField->SetPixel(Index,Shift);
				
				}//End X-Axis
			}//End Y-Axis
		}//End Z-Axis
	}


	//Compute the minimum and maximum intensity for both image..
	template< class TInputImage, class TOutputImage>	
	void Particle< TInputImage, TOutputImage>
	::FindMinMax()
	{

		typename MinMaxFilterType::Pointer  fixedMinMax  = MinMaxFilterType::New();
		typename MinMaxFilterType::Pointer 	movingMinMax = MinMaxFilterType::New();

		//Compute MinMax fixed Image
		fixedMinMax->SetImage( FixedImage );
		fixedMinMax->Compute();
		minimumFixed = fixedMinMax->GetMinimum();		
		maximumFixed = fixedMinMax->GetMaximum();	

		//Compute MinMax moving Image
		movingMinMax->SetImage( MovingImage );
		movingMinMax->Compute();
		minimumMoving = movingMinMax->GetMinimum();		
		maximumMoving = movingMinMax->GetMaximum();	

		minimum =  minimumFixed;
		maximum =  maximumFixed;
		//Compute Minimum and Maximum
		if( minimumMoving < minimum)
			minimum = minimumMoving;
		if( maximumMoving > maximum)
			maximum = maximumMoving;

		ban_minmax = true;
	}

	//Compute the particle Weigth
	template< class TInputImage, class TOutputImage>	
	float Particle< TInputImage, TOutputImage>
	::ComputeParticleWeight()
	{

		if( !ban_minmax )
			FindMinMax();
			
		GenerateTransformationMatrix();

		if( ImageDimension == 3)
			ComputeParticleWeight3D();
		else
			ComputeParticleWeight2D();	

		return Weight;
	}

	//Compute the particle Weigth 2D
	template< class TInputImage, class TOutputImage>	
	void Particle< TInputImage, TOutputImage>
	::ComputeParticleWeight2D()
	{
		float cx,cy;
		float A00,A01,A02;
		float A10,A11,A12;
		float ax,ay;
		float	binsDivMoving;
		float	binsDivFixed;
		int **Hist = NULL;
		int i,cont;
		int indFixed;
		int indMoving;
		IndexType 			Index;
		ContinuousIndex 	NewPosition;
		PixelType 			ValueMoving;	
		PixelType 			ValueFixed;	

		//Compute the size of step according to number of bins for histogram computation
		binsDivFixed = (float)NumberOfBins/(float)(maximumFixed - minimumFixed + 1);
		binsDivMoving = (float)NumberOfBins/(float)(maximumMoving - minimumMoving + 1);

		//Alloacate Memory for Joint Histogram
		Hist = new int *[NumberOfBins];
		for(i = 0; i < NumberOfBins; i++){
			Hist[i] = new int [NumberOfBins];
			memset(Hist[i],0,sizeof(int)*NumberOfBins);
		}

		//Obtain Image Center
		cx = MovingImageSize[0]/2.0;	//X-axis center
		cy = MovingImageSize[1]/2.0;	//Y-axis center

		A02 = AffineMatrix(0,2);	//Dx
		A12 = AffineMatrix(1,2);	//Dy

		//Asign MovingImage to Interpolator
		m_Interpolator->SetInputImage( MovingImage );

		cont = 0;
		//Y-Axis
		for(Index[1] = 0; Index[1] <FixedImageSize[1]; Index[1] += SamplesStepSize[1])
		{
			ay = Index[1] - cy;
			A01 = AffineMatrix(0,1) * ay;
			A11 = AffineMatrix(1,1) * ay;
			//X-Axis
			for(Index[0] = 0; Index[0] < FixedImageSize[0]; Index[0] += SamplesStepSize[0])
			{
				ax = Index[0] - cx;
				A00 = AffineMatrix(0,0) * ax;
				A10 = AffineMatrix(1,0) * ax;

				//Obtain the New Position
				NewPosition[0] = A00 + A01 + A02 + cx;
				NewPosition[1] = A10 + A11 + A12 + cy;

				//Check the Limits of the Image
				if( m_Interpolator->IsInsideBuffer( NewPosition ) )
				{	
					//Use the Interpolator to Obtain MovingImage Value
					ValueMoving = static_cast<PixelType>
							  		( m_Interpolator->EvaluateAtContinuousIndex( NewPosition ) );
//						ValueMoving = static_cast<PixelType>
//								  		( m_Interpolator->Evaluate( Point ) );

					//Get FixedImage Value
					ValueFixed = FixedImage->GetPixel( Index );

					//Compute histogram indexs
					indMoving = (ValueMoving - minimumMoving)*binsDivMoving;	
					indFixed = (ValueFixed - minimumFixed)*binsDivFixed;	
					
					//Construct Joint Hist
					Hist[indMoving][indFixed]++; 
			
					cont++;//Counter of samples
				}
								
			}//End X-Axis
		}//End Y-Axis

		//Compute the similarity 
		SimilarityValue = ComputeSimilarityMetric(Hist,cont);

		c3 = 2.0 - SimilarityValue;

		//Compute Likelihood - Particle Weight
		Weight = c1*exp( c2* c3* c3);

		//Releasing Memory
		for(i = 0; i < NumberOfBins; i++)
			delete [] Hist[i];
		delete [] Hist;

		m_Interpolator->SetInputImage( NULL );

	}
	

	//Compute the particle Weigth 2D
	template< class TInputImage, class TOutputImage>	
	void Particle< TInputImage, TOutputImage>
	::ComputeParticleWeight3D()
	{
		float cx,cy,cz;
		float A00,A01,A02,A03;
		float A10,A11,A12,A13;
		float A20,A21,A22,A23;
		float ax,ay,az;
		float	binsDivMoving;
		float	binsDivFixed;
		int **Hist = NULL;
		int i,cont;
		int indFixed;
		int indMoving;
		IndexType 			Index;
		PointType			Point;
		ContinuousIndex 	Shift;
		ContinuousIndex 	NewPosition;
		PixelType 			ValueMoving;	
		PixelType 			ValueFixed;	

		//Compute the size of step according to number of bins for histogram computation
		binsDivFixed = (float)NumberOfBins/(float)(maximumFixed - minimumFixed + 1);
		binsDivMoving = (float)NumberOfBins/(float)(maximumMoving - minimumMoving + 1);

		//Alloacate Memory for Joint Histogram
		Hist = new int *[NumberOfBins];
		for(i = 0; i < NumberOfBins; i++){
			Hist[i] = new int [NumberOfBins];
			memset(Hist[i],0,sizeof(int)*NumberOfBins);
		}

		//Obtain Image Center
		cx = MovingImageSize[0]/2.0;	//X-axis center
		cy = MovingImageSize[1]/2.0;	//Y-axis center
		cz = MovingImageSize[2]/2.0;	//Z-axis center

		A03 = AffineMatrix(0,3);	//Dx
		A13 = AffineMatrix(1,3);	//Dy
		A23 = AffineMatrix(2,3);	//Dz
		
		//Asign MovingImage to Interpolator
		m_Interpolator->SetInputImage( MovingImage );

		cont = 0;
		//Z-Axis
		for(Index[2] = 0; Index[2] < FixedImageSize[2]; Index[2] += SamplesStepSize[2])
		{
			az = Index[2] - cz;
			A02 = AffineMatrix(0,2) * az;
			A12 = AffineMatrix(1,2) * az;
			A22 = AffineMatrix(2,2) * az;
			//Y-Axis
			for(Index[1] = 0; Index[1] < FixedImageSize[1]; Index[1] += SamplesStepSize[1])
			{
				ay = Index[1] - cy;
				A01 = AffineMatrix(0,1) * ay;
				A11 = AffineMatrix(1,1) * ay;
				A21 = AffineMatrix(2,1) * ay;
				//X-Axis
				for(Index[0] = 0; Index[0] < FixedImageSize[0]; Index[0] += SamplesStepSize[0])
				{
					ax = Index[0] - cx;
					A00 = AffineMatrix(0,0) * ax;
					A10 = AffineMatrix(1,0) * ax;
					A20 = AffineMatrix(2,0) * ax;

					Shift[0] = A00 + A01 + A02 + A03 + cx - Index[0];
					Shift[1] = A10 + A11 + A12 + A13 + cy - Index[1];
					Shift[2] = A20 + A21 + A22 + A23 + cz - Index[2];

					//--------- Physical Coordinates --------------//
					if( m_PhysicalCoordinates )
					{
						MovingImage->TransformIndexToPhysicalPoint(Index,Point);

						for( i = 0; i < ImageDimension; i++)
							Point[i] += Shift[i];

						//Check the Limits of the Image
						if( m_Interpolator->IsInsideBuffer( Point ) )
						{	
							ValueMoving = static_cast<PixelType>
									  		( m_Interpolator->Evaluate( Point ) );
							//Get FixedImage Value
							ValueFixed = FixedImage->GetPixel( Index );

							//Compute histogram indexs
							indMoving = (ValueMoving - minimumMoving)*binsDivMoving;	
							indFixed = (ValueFixed - minimumFixed)*binsDivFixed;	
					
							//Construct Joint Hist
							Hist[indMoving][indFixed]++; 
			
							cont++;//Counter of samples
						}
					}
					//--------- End Physical Coordinates --------------//

					//--------- Voxel Coordinates ---------------//
					else
					{
						for( i = 0; i < ImageDimension; i++)
							NewPosition[i] = Index[i] + Shift[i];

						if( m_Interpolator->IsInsideBuffer( NewPosition ) )
						{	

							//Use the Interpolator to Obtain MovingImage Value
							ValueMoving = static_cast<PixelType>
									  		( m_Interpolator->EvaluateAtContinuousIndex( NewPosition ) );
							//Get FixedImage Value
							ValueFixed = FixedImage->GetPixel( Index );

							//Compute histogram indexs
							indMoving = (ValueMoving - minimumMoving)*binsDivMoving;	
							indFixed = (ValueFixed - minimumFixed)*binsDivFixed;	
					
							//Construct Joint Hist
							Hist[indMoving][indFixed]++; 
			
							cont++;//Counter of samples
						}
					}
					//--------- End Voxel Coordinates --------------//


				}//End X-Axis
			}//End Y-Axis
		}//End Z-Axis

		m_Interpolator->SetInputImage( NULL );

		//Compute the similarity 
		SimilarityValue = ComputeSimilarityMetric(Hist,cont);

		c3 = 2.0 - SimilarityValue;

		//Compute Likelihood - Particle Weight
		Weight = c1*exp( c2* c3* c3);

		//Releasing Memory
		for(i = 0; i < NumberOfBins; i++)
			delete [] Hist[i];
		delete [] Hist;
	}

	//Compute the Similarity Metric (Normalized Mutual Information)
	template< class TInputImage, class TOutputImage>	
	float Particle< TInputImage, TOutputImage>
	::ComputeSimilarityMetric(int **Hist,int samples)
	{
		float jointEntropy,fixedEntropy,movingEntropy;
		float probability,sum,SM;
		int i,j,S1,S2;

		sum = 1.0/(float)samples;

		//Entropies
		jointEntropy = 0.0;
		fixedEntropy = 0.0;
		movingEntropy = 0.0;

		//Compute Entropies
		for(i = 0; i < NumberOfBins; i++)
		{
			S1 = 0;
			S2 = 0;
			for(j = 0; j < NumberOfBins; j++)
			{
				S1 += Hist[i][j];
				S2 += Hist[j][i];
				//Joint Entropy
				if( Hist[i][j] > 0 )
				{
				   probability = (float) Hist[i][j] * sum;
				   jointEntropy -= probability * log2( probability );
				}
			}
			//Moving Entropy
			if( S1 > 0 )
			{
			   probability = (float) S1 * sum;
			   movingEntropy -= probability * log2( probability );
			}
			//Fixed Entropy
			if( S2 > 0 )
			{
			   probability = (float) S2 * sum;
			   fixedEntropy -= probability * log2( probability );
			}		
		}//End first For

		
		SM = (fixedEntropy + movingEntropy) / jointEntropy;
		
		return SM;
	}

	//Get the transformed Image
	template< class TInputImage, class TOutputImage>	
	typename TOutputImage::Pointer Particle< TInputImage, TOutputImage>
	::GetTransformedImage()
	{
		//For class Particle it is not necessary to make something in this function
		typename OutputCasterType::Pointer ImageCaster = OutputCasterType::New();
		typename InternalDuplicatorType::Pointer		Duplicator = InternalDuplicatorType::New();

		//If Vector Field is not sended
		if( !ban_send_vec )
		{
			GenerateDisplacementField();
		}

		FindMinMax();

		//--------------- Voxel Coordinates -------------///
		if( !m_PhysicalCoordinates )
		{
			//Allocate Output Image
			Duplicator->SetInputImage( FixedImage );
			Duplicator->Update();
			AuxMovingImage = Duplicator->GetOutput();
			AuxMovingImage->FillBuffer( minimumMoving );

			//Declarate Iterators		
			itk::ImageRegionIteratorWithIndex<VectorImageType>
							fieldIte( DisplacementField,DisplacementField->GetRequestedRegion() );
			itk::ImageRegionIterator< InternalImageType > 
							outputIte( AuxMovingImage,AuxMovingImage->GetRequestedRegion() );

			//Auxiliar variables
			IndexType 			Index;
			PixelType 			Value;	
			VectorPixelType 	Displacement;
			ContinuousIndex 	NewPosition;

			//Interpolator
			m_Interpolator->SetInputImage( MovingImage );

			//Iterators
			fieldIte.GoToBegin();
			outputIte.GoToBegin();
		
			//Obtain the Transformed Image
			while( !outputIte.IsAtEnd() )
			{
				//GetOutputImageIndex
				Index = outputIte.GetIndex();

				//Get Displacement
				Displacement = fieldIte.Get();	

				//Obtain the new voxel position
				for( int i = 0; i < ImageDimension; i++)
					NewPosition[i] = Index[i] + Displacement[i];
			
				//Check the Limits of the Image
				if( m_Interpolator->IsInsideBuffer( NewPosition ) )
				{	//Use the Interpolator
					Value = static_cast<PixelType>
							  ( m_Interpolator->EvaluateAtContinuousIndex( NewPosition ) );
				}
				else
				{
					Value = minimumMoving;
				}
				
				outputIte.Set( Value );
				++fieldIte;
				++outputIte;
			}

			ImageCaster->SetInput( AuxMovingImage );
		}
		//------------------ End Voxel Coordinates ----------------//

		
		//-------------- Physical Coordinates ---------------//
		else
		{

			typename WarpFilterType::Pointer WarpFilter = WarpFilterType::New();

			//Set Warp Filter Parameters
			WarpFilter->SetDisplacementField( DisplacementField );
			WarpFilter->SetInterpolator( m_Interpolator );
			WarpFilter->SetInput( MovingImage );
			WarpFilter->SetOutputOrigin( FixedImage->GetOrigin() );
			WarpFilter->SetOutputSpacing( FixedImage->GetSpacing() );
			WarpFilter->SetOutputDirection( FixedImage->GetDirection() );
			WarpFilter->SetEdgePaddingValue( minimumMoving );
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

			//Cast InternalImageType to OutputType
			ImageCaster->SetInput( WarpFilter->GetOutput() );
		}

		//------------------ End Physical Coordinates ----------------//

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

		//Return Output Image
		OutputImage = ImageCaster->GetOutput();

		return OutputImage;
	}

	//********************************************
	//Generate Data   (Main fucntion of the Class)
	//********************************************
	template< class TInputImage, class TOutputImage>	
	void Particle< TInputImage, TOutputImage>
	::GenerateData()
	{
		
	}	


}//namespace itk 

#endif
