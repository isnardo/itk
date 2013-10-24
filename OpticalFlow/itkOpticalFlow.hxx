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

		if(ImageDimension == 2){
			m_LambdaXY = 6000;
			m_LambdaZ  = 15000;
			m_LambdaW  = 20000;
			m_Iterations = 5;
			m_GaussSeidelIterations = 200;
		}
		else
		{
			m_LambdaXY = 6000;
			m_LambdaZ  = 15000;
			m_LambdaW  = 20000;
			m_Iterations = 5;
			m_GaussSeidelIterations = 200;			
		}
		m_Threshold = 1.0;

		//Default Interpolator
		typename DefaultInterpolatorType::Pointer Interp = DefaultInterpolatorType::New();
	  m_Interpolator =  static_cast<InterpolatorType*>( Interp.GetPointer() );
	}


	//********************************************
	//Generate Data   (Main fucntion of the Class)
	//********************************************
	template< class TInputImage, class TOutputImage>	
	void OpticalFlow< TInputImage, TOutputImage>
	::GenerateData()
	{

	}	

}//namespace itk 

#endif
