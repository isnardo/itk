#ifndef _itkMultiScaleOpticalFlow_h
#define _itkMultiScaleOpticalFlow_h

#include "itkOpticalFlow.h"

namespace itk
{
        //Template Class Multi Scale Optical Flow
        template< class TInputImage, class TOutputImage>
        class ITK_EXPORT MultiScaleOpticalFlow : public OpticalFlow< TInputImage, TOutputImage >
        {
                public:
                        //** Standard class typedefs **//
                        typedef MultiScaleOpticalFlow                                 Self;
                        typedef OpticalFlow< TInputImage,TOutputImage >        				Superclass;
                        typedef SmartPointer< Self >                                  Pointer;
                        typedef SmartPointer< const Self >                            ConstPointer;

                        itkNewMacro( Self );
                        itkTypeMacro( OpticalFlow,ImageToImageFilter );

                protected:

                        MultiScaleOpticalFlow();
                        ~MultiScaleOpticalFlow(){};

                        virtual void GenerateData();

                private:
                        //Private Variables
 

        };//Class Multi Scale Optical Flow

        
} //namespace ITK

#ifndef ITK_MANUAL_INSTANTIATION
#include "itkMultiScaleOpticalFlow.hxx"
#endif

#endif
