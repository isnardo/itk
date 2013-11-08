#ifndef _itkOpticalFlow_txx
#define _itkOpticalFlow_txx

#include "itkOpticalFlow.h"

namespace itk
{
	//*******************************************
	// Functions definition of Optical Flow Class
	//*******************************************
		
	//Class Builder
	template< class TInputImage, class TOutputImage>	
	OpticalFlow< TInputImage, TOutputImage>
	::OpticalFlow()
	{
		//Number of requaired images to Update
		this->SetNumberOfRequiredInputs( 0 );

		//Default variables values

		if( ImageDimension == 2 ){
			m_LambdaXY = 3000;
			m_LambdaZ  = 0;
			m_LambdaW  = 0;
		}
		else
		{
			m_LambdaXY = 6000;
			m_LambdaZ  = 15000;
			m_LambdaW  = 20000;			
		}

		m_Iterations = 5;
		m_GaussSeidelIterations = 200;
		m_Threshold = 1.0;
		m_ShowIterations = false;
		DisplacementFieldFlag = false;

		//Default Interpolator
		typename DefaultInterpolatorType::Pointer Interp = DefaultInterpolatorType::New();
	  m_Interpolator =  static_cast<InterpolatorType*>( Interp.GetPointer() );
	}

	//Set Fixed Image
	template< class TInputImage, class TOutputImage>	
	void OpticalFlow< TInputImage, TOutputImage>
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
	}

	//Set Moving Image
	template< class TInputImage, class TOutputImage>	
	void OpticalFlow< TInputImage, TOutputImage>
	::SetMovingImage(const TInputImage* movingImage)
	{
		//this->SetNthInput(1, const_cast<TInputImage*>(movingImage));
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
	}

	//Compute Intensities Range Between Images
	template< class TInputImage, class TOutputImage>	
	float OpticalFlow< TInputImage, TOutputImage>
	::ComputeIntensitiesRange()
	{
		typename ImageMinMaxFilterType::Pointer  	fixedMinMax  = ImageMinMaxFilterType::New();
		typename ImageMinMaxFilterType::Pointer 	movingMinMax = ImageMinMaxFilterType::New();

		float minimumFixed;
		float maximumFixed;
		float minimumMoving;
		float maximumMoving;
		float minimum;
		float maximum;

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

		return (maximum-minimum);
	}

	//Compute Moving Image Gradient
	template< class TInputImage, class TOutputImage>	
	void OpticalFlow< TInputImage, TOutputImage>
	::ComputeGradient()
	{
		//Derivadas
		typename ConvolutionFilterType::Pointer		gxFilter = ConvolutionFilterType::New();
		typename ConvolutionFilterType::Pointer		gyFilter = ConvolutionFilterType::New();
		typename ConvolutionFilterType::Pointer		gzFilter = ConvolutionFilterType::New();
		typename InternalImageType::Pointer				Kernel_x = InternalImageType::New();
		typename InternalImageType::Pointer				Kernel_y = InternalImageType::New();
		typename InternalImageType::Pointer				Kernel_z = InternalImageType::New();
		SizeType  																kernel_size;
		IndexType 																kernel_start;
		IndexType 																kernel_index;
		RegionType   															kernel_region;
		PixelType																	kernel_pixel;

		//Define Kernel Region
		kernel_start.Fill(0);
		kernel_size.Fill(3); 
		kernel_pixel = 0.0;
		kernel_region.SetSize( kernel_size );
		kernel_region.SetIndex( kernel_start );	

		//Allocate Kernels Memory
		Kernel_x->SetRegions( kernel_region );
		Kernel_x->Allocate();
		Kernel_x->FillBuffer( kernel_pixel );

		Kernel_y->SetRegions( kernel_region );
		Kernel_y->Allocate();
		Kernel_y->FillBuffer( kernel_pixel );

		Kernel_z->SetRegions( kernel_region );
		Kernel_z->Allocate();
		Kernel_z->FillBuffer( kernel_pixel );

		//2D
		if( ImageDimension == 2)
		{
			//Kernel to derivative X-Axis
			kernel_index[0] = 0;
			kernel_index[1] = 1;
			kernel_pixel = 0.5;
			Kernel_x->SetPixel( kernel_index,kernel_pixel );

			kernel_index[0] = 2;
			kernel_index[1] = 1;
			kernel_pixel = -0.5;
			Kernel_x->SetPixel( kernel_index,kernel_pixel );

			//Kernel to derivative Y-Axis
			kernel_index[0] = 1;
			kernel_index[1] = 0;
			kernel_pixel = 0.5;
			Kernel_y->SetPixel( kernel_index,kernel_pixel );

			kernel_index[0] = 1;
			kernel_index[1] = 2;
			kernel_pixel = -0.5;
			Kernel_y->SetPixel( kernel_index,kernel_pixel );
		}
		//Volume
		else	
		{
			//Kernel to derivative X-Axis
			kernel_index[0] = 0;
			kernel_index[1] = 1;
			kernel_index[2] = 1;
			kernel_pixel = 0.5;
			Kernel_x->SetPixel( kernel_index,kernel_pixel );

			kernel_index[0] = 2;
			kernel_index[1] = 1;
			kernel_index[2] = 1;
			kernel_pixel = -0.5;
			Kernel_x->SetPixel( kernel_index,kernel_pixel );

			//Kernel to derivative Y-Axis
			kernel_index[0] = 1;
			kernel_index[1] = 0;
			kernel_index[2] = 1;
			kernel_pixel = 0.5;
			Kernel_y->SetPixel( kernel_index,kernel_pixel );

			kernel_index[0] = 1;
			kernel_index[1] = 2;
			kernel_index[2] = 1;
			kernel_pixel = -0.5;
			Kernel_y->SetPixel( kernel_index,kernel_pixel );

			//Kernel to derivate Z-Axis
			kernel_index[0] = 1;
			kernel_index[1] = 1;
			kernel_index[2] = 0;
			kernel_pixel = 0.5;
			Kernel_z->SetPixel( kernel_index,kernel_pixel );

			kernel_index[0] = 1;
			kernel_index[1] = 1;
			kernel_index[2] = 2;
			kernel_pixel = -0.5;
			Kernel_z->SetPixel( kernel_index,kernel_pixel );

		}

		//Derivative in X-Axis
		gxFilter->SetInput( MovingImage );
		gxFilter->SetKernelImage( Kernel_x );
		try
		{
			gxFilter->Update();
		}
		catch( itk::ExceptionObject & excep )
		{
			std::cerr << "Exception catched !" << std::endl;
			std::cerr << excep << std::endl;
		}
		DerivativeX = gxFilter->GetOutput();

		//Derivative in Y-Axis
		gyFilter->SetInput( MovingImage );
		gyFilter->SetKernelImage( Kernel_y);
		try
		{
			gyFilter->Update();
		}
		catch( itk::ExceptionObject & excep )
		{
			std::cerr << "Exception catched !" << std::endl;
			std::cerr << excep << std::endl;
		}
		DerivativeY = gyFilter->GetOutput();

		//Volume
		if( ImageDimension == 3){
			//Derivative in Z-Axis
			gzFilter->SetInput( MovingImage );
			gzFilter->SetKernelImage( Kernel_z);
			try
			{
				gzFilter->Update();
			}
			catch( itk::ExceptionObject & excep )
			{
				std::cerr << "Exception catched !" << std::endl;
				std::cerr << excep << std::endl;
			}
			DerivativeZ = gzFilter->GetOutput();
		}	
	}

	//Allocate Memory for Displacement Field
	template< class TInputImage, class TOutputImage>	
	void OpticalFlow< TInputImage, TOutputImage>
	::AllocateDisplacementField()
	{
		VectorPixelType 	valueFill;

		DisplacementField = VectorImageType::New();
		DisplacementField->SetSpacing( FixedImage->GetSpacing() );
		DisplacementField->SetOrigin( FixedImage->GetOrigin() );
		DisplacementField->SetDirection( FixedImage->GetDirection() );
		DisplacementField->SetRegions( FixedImage->GetLargestPossibleRegion() );
		DisplacementField->Allocate();		

		valueFill.Fill( 0.0 );
		DisplacementField->FillBuffer( valueFill );
	}

	//Compute Percentage Differences Between the Images Intensities
	template< class TInputImage, class TOutputImage>	
	float OpticalFlow< TInputImage, TOutputImage>
	::ComputePercentageImageDifferences()
	{
		itk::ImageRegionIterator<InternalImageType>  it( SubImage,SubImage->GetRequestedRegion() );

		float 	s = 0.0;
		int 		k = 0;

		IndexType Index;
		it.GoToBegin();

		if( ImageDimension == 3 )//3D
			//Differences in Intensisties
			while( !it.IsAtEnd() )
			{
				Index = it.GetIndex();
				if( Index[2] > 0 && Index[2] < FixedImageSize[2]-1 )
				{
					s += sqrt( it.Get()*it.Get() );
					k++; 
				}
				++it;
			}		
		else//2D
			//Differences in Intensisties
			while( !it.IsAtEnd() )
			{
				s += sqrt( it.Get()*it.Get() );
				k++; 
				++it;
			}		

		s /= (float)k;
		return s*100.0/this->IntensitiesRange;
	}

	//Compute Optical Flow 
	template< class TInputImage, class TOutputImage>	
	void OpticalFlow< TInputImage, TOutputImage>
	::ComputeOpticalFlow()
	{
		typename AddFilterType::Pointer 				AddFilter = AddFilterType::New();
		typename SubFilterType::Pointer 				SubFilter = SubFilterType::New();


		//Apply Vector Field to Derivatives
		aux_DerivativeX = ApplyDisplacementField( DerivativeX,DisplacementField );
		aux_DerivativeY = ApplyDisplacementField( DerivativeY,DisplacementField );
		//Volume ( 3D )
		if( ImageDimension == 3)
		{
			aux_DerivativeZ = ApplyDisplacementField( DerivativeZ,DisplacementField );
		}
	
		SubFilter->SetInput1( FixedImage );
		SubFilter->SetInput2( ApplyDisplacementField( MovingImage,DisplacementField ) );
		try
		{
			SubFilter->Update();
		}
		catch( itk::ExceptionObject & excep )
		{
			std::cerr << "Exception catched !" << std::endl;
			std::cerr << excep << std::endl;
		}
		SubImage = SubFilter->GetOutput();
		
		//Volume (3D)
		if( ImageDimension == 3 )
		{
			//Start Gauss Seidel Method to Solve the Linear Equations System for a Volume (3D images)
			this->SolveGaussSeidel3D();
		}
		//2D
		else
		{	
			//Start Gauss Seidel Method to Solve the Linear Equations System for 2D images
			this->SolveGaussSeidel2D();
		}

		//Add Rigid Vector Field to Optical Flow Vector Field
		AddFilter->SetInput1( DisplacementField );
		AddFilter->SetInput2( aux_DisplacementField );
		try
		{
			AddFilter->Update();
		}
		catch( itk::ExceptionObject & excep )
		{
			std::cerr << "Exception catched !" << std::endl;
			std::cerr << excep << std::endl;
		}

		//Obtain the Non-Rigid Deformation Vector Field
		DisplacementField = AddFilter->GetOutput();

		//Update Differences between Images
		SubFilter->SetInput1( FixedImage );
		SubFilter->SetInput2( ApplyDisplacementField( MovingImage,DisplacementField ) );
		try
		{
			SubFilter->Update();
		}
		catch( itk::ExceptionObject & excep )
		{
			std::cerr << "Exception catched !" << std::endl;
			std::cerr << excep << std::endl;
		}
		SubImage = SubFilter->GetOutput();
	}


	//Function to return moving image with displacement field applied
	template< class TInputImage, class TOutputImage>	
	typename TOutputImage::Pointer OpticalFlow< TInputImage, TOutputImage>
	::GetTransformedImage()
	{
		typename OutputCasterType::Pointer  ImageCaster = OutputCasterType::New();

		//Apply Vdisplacement Field
		ImageCaster->SetInput( ApplyDisplacementField( MovingImage,DisplacementField ) );

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


	//Compute Gauss Seidel to Solve Linear System of Equations (2D)
	template< class TInputImage, class TOutputImage>	
	void OpticalFlow< TInputImage, TOutputImage>
	::SolveGaussSeidel2D()
	{
		//Iterators
		itk::ImageRegionIteratorWithIndex< VectorImageType > 
														IteVec( aux_DisplacementField,aux_DisplacementField->GetRequestedRegion() );
		itk::ImageRegionIterator< InternalImageType > 
														IteGx( aux_DerivativeX,aux_DerivativeX->GetRequestedRegion() );
		itk::ImageRegionIterator< InternalImageType > 
														IteGy( aux_DerivativeY,aux_DerivativeY->GetRequestedRegion() );
		itk::ImageRegionIterator< InternalImageType > 
														IteSub( SubImage,SubImage->GetRequestedRegion() );

		IndexType								Index;
		PixelType								Value;	
		IndexType								IndexG;
		PixelType								ValueGx;
		PixelType								ValueGy;
		VectorIndexType 				IndexVec;
		VectorPixelType 				ValueVec;
		VectorPixelType 				ValueVecA;
		VectorPixelType 				valueFill;

		float Sx;
		float Sy;
		int   Ns;

		//Set Zero Vector Field
		valueFill.Fill( 0.0 );
		aux_DisplacementField->FillBuffer( valueFill );

		//Gauss-Seidel
		for(int i = 0; i < m_GaussSeidelIterations; i++)
		{
			IteSub.GoToBegin();
			IteGx.GoToBegin();
			IteGy.GoToBegin();
			IteVec.GoToBegin();
			while( !IteVec.IsAtEnd() )
			{
				Sx = 0.0;
				Sy = 0.0;
				Ns = 0;		//Cardinality
			
				//Z-Axis Neihgborhood
				Index = IteVec.GetIndex();
				IndexVec[0] = Index[0];
				if(Index[1] > 0)// (x,y-1)
				{	
					IndexVec[1] = Index[1] - 1;
					ValueVec = aux_DisplacementField->GetPixel( IndexVec );
					Sx += ValueVec[0];
					Sy += ValueVec[1];
					Ns++;
				}
				if(Index[1] < FixedImageSize[1]-1)// (x,y+1)
				{
					IndexVec[1] = Index[1] + 1;
					ValueVec = aux_DisplacementField->GetPixel( IndexVec );
					Sx += ValueVec[0];
					Sy += ValueVec[1];
					Ns++;
				}
				//X-Axis Neihgborhood
				IndexVec[1] = Index[1];
				if(Index[0] > 0)// (x-1,y)
				{	
					IndexVec[0] = Index[0] - 1;
					ValueVec = aux_DisplacementField->GetPixel( IndexVec );
					Sx += ValueVec[0];
					Sy += ValueVec[1];
					Ns++;
				}
				if(Index[0] < FixedImageSize[0]-1)// (x+1,y)
				{
					IndexVec[0] = Index[0] + 1;
					ValueVec = aux_DisplacementField->GetPixel( IndexVec );
					Sx += ValueVec[0];
					Sy += ValueVec[1];
					Ns++;
				}

				Value = IteSub.Get();
				ValueGx = IteGx.Get();
				ValueGy = IteGy.Get();
				ValueVecA = IteVec.Get();

				//Linear System of Equations

				//Displacement X-axis
				ValueVec[0] = ( ValueGx*( Value - ValueGy*ValueVecA[1] ) + m_LambdaXY*Sx )
									/( ValueGx*ValueGx + m_LambdaXY*Ns + m_LambdaW );
				//Displacement Y-axis
				ValueVec[1] = ( ValueGy*( Value - ValueGx*ValueVecA[0] ) + m_LambdaXY*Sy )
									/( ValueGy*ValueGy + m_LambdaXY*Ns + m_LambdaW );

				IteVec.Set( ValueVec );

				++IteSub;
				++IteGx;
				++IteGy;
				++IteVec;
			}//End While
		}// End Gauss-Seidel

	}


	//Compute Gauss Seidel to Solve Linear System of Equations (3D)
	template< class TInputImage, class TOutputImage>	
	void OpticalFlow< TInputImage, TOutputImage>
	::SolveGaussSeidel3D()
	{

		//Iterators
		itk::ImageRegionIteratorWithIndex< VectorImageType > 
														IteVec( aux_DisplacementField,aux_DisplacementField->GetRequestedRegion() );
		itk::ImageRegionIterator< InternalImageType > 
														IteGx( aux_DerivativeX,aux_DerivativeX->GetRequestedRegion() );
		itk::ImageRegionIterator< InternalImageType > 
														IteGy( aux_DerivativeY,aux_DerivativeY->GetRequestedRegion() );
		itk::ImageRegionIterator< InternalImageType > 
														IteGz( aux_DerivativeZ,aux_DerivativeZ->GetRequestedRegion() );
		itk::ImageRegionIterator< InternalImageType > 
														IteSub( SubImage,SubImage->GetRequestedRegion() );

		//Local Variables
		IndexType 								Index;
		PixelType 								Value;	
		PixelType									ValueGx;
		PixelType									ValueGy;
		PixelType									ValueGz;
		VectorIndexType 					IndexVec;
		VectorPixelType 					ValueVec;
		VectorPixelType 					ValueVecA;
		VectorPixelType 					valueFill;

		float 		Sx;
		float 		Sy;
		float 		Sz;
		int   		Ns;

		//Set Zero Vector Field
		valueFill.Fill( 0.0 );
		aux_DisplacementField->FillBuffer( valueFill );

		//Gauss-Seidel
		for(int i = 0; i < m_GaussSeidelIterations; i++)
		{
			IteSub.GoToBegin();
			IteGx.GoToBegin();
			IteGy.GoToBegin();
			IteGz.GoToBegin();
			IteVec.GoToBegin();
			while( !IteVec.IsAtEnd() )
			{
				Sx = 0.0;
				Sy = 0.0;
				Sz = 0.0;
				Ns = 0;		//Cardinality
			
				//Z-Axis Neihgborhood
				Index = IteVec.GetIndex();
				IndexVec[0] = Index[0];
				IndexVec[1] = Index[1];
				if(Index[2] > 0)// (x,y,z-1)
				{
					IndexVec[2] = Index[2] - 1;
					ValueVec = aux_DisplacementField->GetPixel( IndexVec );
					Sx += ValueVec[0];
					Sy += ValueVec[1];
					Sz += ValueVec[2];
					Ns++;
				}
				if(Index[2] < FixedImageSize[2]-1)// (x,y,z+1)
				{
					IndexVec[2] = Index[2] + 1;
					ValueVec = aux_DisplacementField->GetPixel( IndexVec );
					Sx += ValueVec[0];
					Sy += ValueVec[1];
					Sz += ValueVec[2];
					Ns++;
				}
				//Y-Axis Neihgborhood
				IndexVec[0] = Index[0];
				IndexVec[2] = Index[2];
				if(Index[1] > 0)// (x,y-1,z)
				{
					IndexVec[1] = Index[1] - 1;
					ValueVec = aux_DisplacementField->GetPixel( IndexVec );
					Sx += ValueVec[0];
					Sy += ValueVec[1];
					Sz += ValueVec[2];
					Ns++;
				}
				if(Index[1] < FixedImageSize[1]-1){//x,y+1
					IndexVec[1] = Index[1] + 1;
					ValueVec = aux_DisplacementField->GetPixel( IndexVec );
					Sx += ValueVec[0];
					Sy += ValueVec[1];
					Sz += ValueVec[2];
					Ns++;
				}
				//X-Axis Neihgborhood
				IndexVec[1] = Index[1];
				IndexVec[2] = Index[2];
				if(Index[0] > 0)// (x-1,y,z)
				{	
					IndexVec[0] = Index[0] - 1;
					ValueVec = aux_DisplacementField->GetPixel( IndexVec );
					Sx += ValueVec[0];
					Sy += ValueVec[1];
					Sz += ValueVec[2];
					Ns++;
				}
				if(Index[0] < FixedImageSize[0]-1)// (x+1,y,z)
				{
					IndexVec[0] = Index[0] + 1;
					ValueVec = aux_DisplacementField->GetPixel( IndexVec );
					Sx += ValueVec[0];
					Sy += ValueVec[1];
					Sz += ValueVec[2];
					Ns++;
				}

				Value = IteSub.Get();
				ValueGx = IteGx.Get();
				ValueGy = IteGy.Get();
				ValueGz = IteGz.Get();
				ValueVecA = IteVec.Get();

				//Linear System of Equations

				//Displacement X-axis
				ValueVec[0] = ( ValueGx*( Value - ValueGy*ValueVecA[1] - ValueGz*ValueVecA[2]) 
									+ m_LambdaXY*Sx )
									/( ValueGx*ValueGx + m_LambdaXY*Ns + m_LambdaW );
				//Displacement Y-axis
				ValueVec[1] = ( ValueGy*( Value - ValueGx*ValueVecA[0] - ValueGz*ValueVecA[2] ) 
									+ m_LambdaXY*Sy )
									/( ValueGy*ValueGy + m_LambdaXY*Ns + m_LambdaW );
				//Displacement Z-axis
				ValueVec[2] = ( ValueGz*( Value - ValueGx*ValueVecA[0] - ValueGy*ValueVecA[1]) 
									+ m_LambdaZ*Sz )
									/( ValueGz*ValueGz + m_LambdaZ*Ns + m_LambdaW );

				//OF_VectorField->SetPixel( IndexVec, ValueVec );
				IteVec.Set( ValueVec );

				++IteSub;
				++IteGx;
				++IteGy;
				++IteGz;
				++IteVec;
			}//End While
		}// End Gauss-Seidel
	}

	//********************************************
	//Generate Data   (Main fucntion of the Class)
	//********************************************
	template< class TInputImage, class TOutputImage>	
	void OpticalFlow< TInputImage, TOutputImage>
	::GenerateData()
	{

		typename DuplicatorType::Pointer 	duplicator = DuplicatorType::New();
		float 														error, error_1, dif;

		//Allocate Displacement Field Memory
		if( !DisplacementFieldFlag )
			this->AllocateDisplacementField();

		duplicator->SetInputImage( DisplacementField );
		try
		{
			duplicator->Update();
		}
		catch( itk::ExceptionObject & excep )
		{
			std::cerr << "Exception catched !" << std::endl;
			std::cerr << excep << std::endl;
		}
		aux_DisplacementField = duplicator->GetOutput();
		MovingImage->SetSpacing( FixedImage->GetSpacing() );

		//Compute the Gradient
		this->ComputeGradient();

		//Obtain the Intensities Range
		this->IntensitiesRange = this->ComputeIntensitiesRange();

		//Iterative Optical FLow
		for(int i = 0; i < m_Iterations; i++)
		{

			error_1 = this->ComputePercentageImageDifferences();		

			//Compute The proposed Optical Flow
			this->ComputeOpticalFlow();
/*
///******************
			if( i == 0)
			{
				VectorPixelType 					Value;
				itk::ImageRegionIteratorWithIndex< VectorImageType > IteVec( DisplacementField,DisplacementField->GetRequestedRegion() );
				typename InternalImageType::SpacingType	Space;
				Space = FixedImage->GetSpacing();

				IteVec.GoToBegin();
				while( !IteVec.IsAtEnd() )
				{
					Value = IteVec.Get();
					for(int j = 0; j < ImageDimension; j++)
						Value[j] *= Space[j];
					IteVec.Set( Value );
					++IteVec;
				}
			}
///*******************/
			
			//Compute the percentage error in the intensities
			error = this->ComputePercentageImageDifferences();

			dif = error_1 - error;

			
			
			//Check threshold
			if( error < m_Threshold )
			{
				if( m_ShowIterations )
					std::cout<<" Iteration = "<< i+1 <<" %Diff. = "<<error_1<<";";
				break;
			}

			if ( dif < 0)
			{
				if( m_ShowIterations )
				{
					if (i == 0)
						std::cout<<" Iteration = "<< i+1 <<" %Diff. = "<<error_1<<";";
					else
						std::cout<<" Iteration = "<< i <<" %Diff. = "<<error_1<<";";
				}
				break;
			}

		}

	}	

}//namespace itk 

#endif
