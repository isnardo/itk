#ifndef _itkMultiScaleOpticalFlow_txx
#define _itkMultiScaleOpticalFlow_txx

#include "itkMultiScaleOpticalFlow.h"

namespace itk
{
	//*******************************************
	// Functions definition of MultiScale Optical Flow Class
	//*******************************************
		
	//Class Builder
	template< class TInputImage, class TOutputImage>	
	MultiScaleOpticalFlow< TInputImage, TOutputImage>
	::MultiScaleOpticalFlow()
	{
		//Number of requaired images to Update
		this->SetNumberOfRequiredInputs( 0 );

		//Default Interpolator
		typename DefaultInterpolatorType::Pointer Interp = DefaultInterpolatorType::New();
	  m_Interpolator =  static_cast<InterpolatorType*>( Interp.GetPointer() );

		m_NumberOfScales = 4;
		m_StartScale = 0;
		Iterations = new int [ m_NumberOfScales ];

		Iterations[0] = 20;
		Iterations[1] = 20;
		Iterations[2] = 20;
		Iterations[3] = 20;

		m_ShowIterations = false;
		m_ShowLevels = false;


		if( ImageDimension == 2 ){
			this->m_LambdaXY = 1;
			this->m_LambdaZ  = 0;
			this->m_LambdaW  = 0;
		}
		else
		{
			this->m_LambdaXY = 1;
			this->m_LambdaZ  = 1;
			this->m_LambdaW  = 1;			
		}

		this->m_Threshold = 1.0;
		
	}


	//Allocate Memory for Displacements Field
	template< class TInputImage, class TOutputImage>	
	void MultiScaleOpticalFlow< TInputImage, TOutputImage>
	::AllocateDisplacementField(InternalImageType *image)
	{
		VectorPixelType 	valueFill;

		DisplacementField = VectorImageType::New();
		DisplacementField->SetSpacing( image->GetSpacing() );
		DisplacementField->SetOrigin( image->GetOrigin() );
		DisplacementField->SetDirection( image->GetDirection() );
		DisplacementField->SetRegions( image->GetLargestPossibleRegion() );
		DisplacementField->Allocate();		

		valueFill.Fill( 0.0 );
		DisplacementField->FillBuffer( valueFill );
	}

	//Return Deformed Image
	template< class TInputImage, class TOutputImage>	
	typename TOutputImage::Pointer MultiScaleOpticalFlow< TInputImage, TOutputImage>
	::GetTransformedImage()
	{
		typename InputCasterType::Pointer   InputCaster = InputCasterType::New();
		typename OutputCasterType::Pointer  OutputCaster = OutputCasterType::New();

		//Apply Vdisplacement Field
		InputCaster->SetInput( const_cast<InputImageType*> ( this->GetInput(1) ) );

		//Update Input Caster
		try
		{
			InputCaster->Update();
		}
		catch( itk::ExceptionObject & excep )
		{
			std::cerr << "Exception catched !" << std::endl;
			std::cerr << excep << std::endl;
		}

		OutputCaster->SetInput( ApplyDisplacementField( InputCaster->GetOutput(),DisplacementField ) );

		//Update Output Caster
		try
		{
			OutputCaster->Update();
		}
		catch( itk::ExceptionObject & excep )
		{
			std::cerr << "Exception catched !" << std::endl;
			std::cerr << excep << std::endl;
		}

		return OutputCaster->GetOutput();

	}

	//Move to _Next Level in the Pyram
	template< class TInputImage, class TOutputImage>	
	void MultiScaleOpticalFlow< TInputImage, TOutputImage>
	::MoveToNextLevelInThePyramid(InternalImageType *Previous, InternalImageType* NextLevel)
	{
		typename VectorExpandType::Pointer VectorExpand = VectorExpandType::New();
		itk::ImageRegionIterator< VectorImageType >	Iterator( aux_DisplacementField,aux_DisplacementField->GetRequestedRegion() );	
		typename InternalImageType::SpacingType		Space;

		typename VectorImageType::PixelType 		PixelVec;
		typename VectorImageType::PixelType 		PixelVec1;
		typename VectorImageType::SizeType			SizeVec1;
		typename VectorImageType::SizeType 			SizeVec2;
		typename VectorImageType::PixelType 			SizeVec3;

		//Obtain Size of Previous and Next Level
		SizeVec1 = Previous->GetLargestPossibleRegion().GetSize();
		SizeVec2 = NextLevel->GetLargestPossibleRegion().GetSize();

		//Obtain the scale value to move the next level
		for(int j = 0; j < ImageDimension; j++)
			SizeVec3[j] = 	(float)SizeVec2[j] / (float)SizeVec1[j];		

		//Scale the values of the displacements
		Iterator.GoToBegin();
		while( !Iterator.IsAtEnd() )
		{
			PixelVec = Iterator.Get();
			for(int j =  0; j < ImageDimension; j++)
				PixelVec1[j] = PixelVec[j]*SizeVec3[j];
			Iterator.Set( PixelVec1 );	
			++Iterator;				
		}

		//Go to Next Level
		VectorExpand->SetExpandFactors( SizeVec3 );
		VectorExpand->SetInput( aux_DisplacementField );

		try
		{
			VectorExpand->Update();		
		}
		catch( itk::ExceptionObject & excep )
		{
			std::cerr << "Exception catched !" << std::endl;
			std::cerr << excep << std::endl;
		}
		
		//Allocate Displacemnts Field for the Next level
		DisplacementField = VectorImageType::New();
		DisplacementField = VectorExpand->GetOutput();
		DisplacementField->SetSpacing( NextLevel->GetSpacing() );
		DisplacementField->SetOrigin( NextLevel->GetOrigin() );
		DisplacementField->SetDirection( NextLevel->GetDirection() );	
	}

	//Compute Optical Flow and obtain the vector Fiel at Level-i
	template< class TInputImage, class TOutputImage>	
	void MultiScaleOpticalFlow< TInputImage, TOutputImage>
	::ComputeOpticalFlow(InternalImageType *Fixed,InternalImageType *Moving, int level)
	{
		typename OpticalFlowType::Pointer OpticalFlow = OpticalFlowType::New();
		typename InternalImageType::SpacingType	Space;
		itk::ImageRegionIterator< VectorImageType >	Iterator( DisplacementField,DisplacementField->GetRequestedRegion() );

		VectorPixelType 		PixelVec;
		VectorPixelType 		PixelVec1;

		if( level < m_NumberOfScales-1 )
		{
			
			for(int j = 0; j < ImageDimension; j++)
			{
				Space[j] = 1.0;
			}
			Fixed->SetSpacing( Space );
		}
		else
		{
			Space = Fixed->GetSpacing();
			Iterator.GoToBegin();
			while( !Iterator.IsAtEnd() )
			{
				PixelVec = Iterator.Get();
				for(int j =  0; j < ImageDimension; j++)
				{
					PixelVec1[j] = PixelVec[j]*Space[j];
				}
				Iterator.Set( PixelVec1 );	
				++Iterator;				
			}
		}


		//Set Parameters to Optical Flow Filter
		OpticalFlow->SetInternalFixedImage( Fixed);
		OpticalFlow->SetInternalMovingImage( Moving );
		OpticalFlow->SetDisplacementField( DisplacementField );

		OpticalFlow->SetIterations( Iterations[ level ] );
		OpticalFlow->SetInterpolator( m_Interpolator );
		OpticalFlow->SetThreshold( this->m_Threshold );

		OpticalFlow->SetLambdaXY( this->m_LambdaXY );
		OpticalFlow->SetLambdaZ( this->m_LambdaZ );
		OpticalFlow->SetLambdaW( this->m_LambdaW );

		if( m_ShowIterations )
			OpticalFlow->ShowIterationsOn();

		//compute the Optical flow in this Scale
		try
		{
			OpticalFlow->Update();
		}
		catch( itk::ExceptionObject & excep )
		{
			std::cerr << "Exception catched !" << std::endl;
			std::cerr << excep << std::endl;
		}		

		//Obtain Displacemnts Field			
		aux_DisplacementField = VectorImageType::New();
		aux_DisplacementField = OpticalFlow->GetDisplacementField();

	}


	//********************************************
	//Generate Data   (Main fucntion of the Class)
	//********************************************
	template< class TInputImage, class TOutputImage>	
	void MultiScaleOpticalFlow< TInputImage, TOutputImage>
	::GenerateData()
	{
		typename PyramidFilterType::Pointer PyramidFixed 		= PyramidFilterType::New();
		typename PyramidFilterType::Pointer PyramidMoving		= PyramidFilterType::New();
		typename NormalizeFilterType::Pointer 	normalizeFixed = NormalizeFilterType::New();
		typename NormalizeFilterType::Pointer 	normalizeMoving = NormalizeFilterType::New();
		typename InputImageType::Pointer fixed  =  const_cast<InputImageType*> ( this->GetInput(0) );
		typename InputImageType::Pointer moving =  const_cast<InputImageType*> ( this->GetInput(1) );

		normalizeFixed->SetInput( fixed );
		normalizeMoving->SetInput( moving );

		//Update Normalize fixed
		try
		{
			normalizeFixed->Update();		
		}
		catch( itk::ExceptionObject & excep )
		{
			std::cerr << "Exception catched !" << std::endl;
			std::cerr << excep << std::endl;
		}
		//Update Normalize moving
		try
		{
			normalizeMoving->Update();		
		}
		catch( itk::ExceptionObject & excep )
		{
			std::cerr << "Exception catched !" << std::endl;
			std::cerr << excep << std::endl;
		}

		//Construct Pyramidal Levels for Fixed Image
		PyramidFixed->SetNumberOfLevels( m_NumberOfScales );
		PyramidFixed->SetInput( normalizeFixed->GetOutput() );
		try
		{
			PyramidFixed->Update();		
		}
		catch( itk::ExceptionObject & excep )
		{
			std::cerr << "Exception catched !" << std::endl;
			std::cerr << excep << std::endl;
		}

		//Construct Pyramidal Levels for Moving Image
		PyramidMoving->SetNumberOfLevels( m_NumberOfScales );
		PyramidMoving->SetInput( normalizeMoving->GetOutput() );
		try
		{
			PyramidMoving->Update();		
		}
		catch( itk::ExceptionObject & excep )
		{
			std::cerr << "Exception catched !" << std::endl;
			std::cerr << excep << std::endl;
		}

		//Allocate Displacement Field at the first level of the pyramid
		AllocateDisplacementField( PyramidMoving->GetOutput( m_StartScale ) );

		//NumberOfScales
		for(int i = m_StartScale; i < m_NumberOfScales; i++)
		{
			//std::cout<<std::endl<<"/////*************************************/////"<<std::endl;
			if( m_ShowLevels )
				std::cout<<std::endl<<" SCALE: "<<i+1<<"    Size: "<<PyramidFixed->GetOutput(i)->GetLargestPossibleRegion().GetSize()<<std::endl;
		  
			//Compute Optical Flow and obtain the vector Field at Level-i
			ComputeOpticalFlow( PyramidFixed->GetOutput(i),PyramidMoving->GetOutput(i),i );

			//Go to Low Level in the Pyramid
			if( i < m_NumberOfScales-1 )
			{
				MoveToNextLevelInThePyramid( PyramidFixed->GetOutput(i),PyramidFixed->GetOutput(i+1) );
			}	
		} //End Scales
	}	

}//namespace itk 

#endif
