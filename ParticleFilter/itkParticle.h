#ifndef _itkParticle_h
#define _itkParticle_h

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
        //Template Class Particle for Particle Filter
        template< class TInputImage, class TOutputImage>
        class ITK_EXPORT Particle : public ImageToImageFilter< TInputImage, TOutputImage >
        {
                public:
                        //** Standard class typedefs **//
                	typedef Particle							                                Self;
                	typedef ImageToImageFilter< TInputImage,TOutputImage > 				Superclass;
                  typedef SmartPointer< Self >                                  Pointer;
                  typedef SmartPointer< const Self >                            ConstPointer;

                  itkNewMacro( Self );
                  itkTypeMacro( Particle,ImageToImageFilter );

                  //Type to manage an internal image type
									typedef	TInputImage																						InputImageType;
									typedef	TOutputImage																					OutputImageType;

									//Const Macro to Define Images Dimensions
									itkStaticConstMacro(ImageDimension, unsigned int,InputImageType::ImageDimension);

                  typedef itk::Image< float,itkGetStaticConstMacro(ImageDimension)>	   	InternalImageType;
                  typedef LinearInterpolateImageFunction<InternalImageType,double>			DefaultInterpolatorType;
                  typedef InterpolateImageFunction<InternalImageType,double>            InterpolatorType;

                protected:

                        Particle();
                        ~Particle(){};

                        virtual void GenerateData();

                private:
                        //Private Variables

        };//Class Multi Scale Optical Flow

        
} //namespace ITK

#ifndef ITK_MANUAL_INSTANTIATION
#include "itkParticle.hxx"
#endif

#endif
