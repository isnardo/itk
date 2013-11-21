#ifndef _itkMultiScaleOpticalFlow_h
#define _itkMultiScaleOpticalFlow_h

#include "itkOpticalFlow.h"
#include "itkMultiResolutionPyramidImageFilter.h"
#include "itkVectorExpandImageFilter.h"
#include "itkNormalizeImageFilter.h"
#include "itkImageRegionIterator.h"


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

                  //Type to manage an internal image type
									typedef	TInputImage																										InputImageType;
									typedef	TOutputImage																									OutputImageType;

									//Const Macro to Define Images Dimensions
									itkStaticConstMacro(ImageDimension, unsigned int,InputImageType::ImageDimension);

                  typedef itk::Image< float,itkGetStaticConstMacro(ImageDimension)>	   	InternalImageType;
                  typedef LinearInterpolateImageFunction<InternalImageType,double>			DefaultInterpolatorType;
                  typedef InterpolateImageFunction<InternalImageType,double>            InterpolatorType;
									typedef OpticalFlow<InternalImageType,InternalImageType>            	OpticalFlowType;
									typedef itk::MultiResolutionPyramidImageFilter<InternalImageType,InternalImageType>				PyramidFilterType;
									typedef itk::NormalizeImageFilter< TInputImage,InternalImageType >		NormalizeFilterType;

                  typedef itk::Vector< float,itkGetStaticConstMacro(ImageDimension)>    		VectorPixelType;
                  //Vector Image Type for the displacements vector field
                  typedef itk::Image< VectorPixelType,itkGetStaticConstMacro(ImageDimension)> VectorImageType;
									//Expand Filter to Expand Vector Field to the Next Level
									typedef itk::VectorExpandImageFilter<VectorImageType,VectorImageType> 			VectorExpandType;

									typedef itk::CastImageFilter<InternalImageType,OutputImageType>  OutputCasterType;
									typedef itk::CastImageFilter<InputImageType,InternalImageType>  InputCasterType;


									//Macro to Set other type of interpolator
                  itkSetObjectMacro( Interpolator,InterpolatorType );
									//MAcro to number of scales
                  itkSetMacro( NumberOfScales,int );
                  itkGetMacro( NumberOfScales,int );
									//MAcro Start scale
                  itkSetMacro( StartScale,int );
                  itkGetMacro( StartScale,int );
                  //Show algorithm iterations
                  itkSetMacro( ShowIterations,bool );
                  itkBooleanMacro( ShowIterations );
                  //Show algorithm iterations
                  itkSetMacro( ShowLevels,bool );
                  itkBooleanMacro( ShowLevels );

									//Input Images
                  virtual void SetFixedImage ( const InputImageType *fixedImage )
																						{
																								this->SetNthInput(0, const_cast<InputImageType*>(fixedImage));
																						};
                 virtual void SetMovingImage( const InputImageType *movingImage )
																						{
																								this->SetNthInput(1, const_cast<InputImageType*>(movingImage));
																						};

								//Set Iterations for each scale
								virtual void SetIterationsOfEachScale(int *ite)
																						{
																								Iterations = ite;
																						};

								 virtual typename OutputImageType::Pointer GetTransformedImage();

                 virtual typename VectorImageType::Pointer GetDisplacementField()
                                                                 {
                                                                			return DisplacementField;
                                                                 };
                              

                protected:

                        MultiScaleOpticalFlow();
                        ~MultiScaleOpticalFlow(){};

                        virtual void GenerateData();

                private:
                        //Private Variables
												int 			m_NumberOfScales;
												int 			*Iterations;
												int				m_StartScale;
												bool			m_ShowLevels;
												bool			m_ShowIterations;
												



												typename InterpolatorType::Pointer 				m_Interpolator;
												typename VectorImageType::Pointer         DisplacementField;
												typename VectorImageType::Pointer         aux_DisplacementField;

												void AllocateDisplacementField(InternalImageType *image);
												void ComputeOpticalFlow(InternalImageType *fixed,InternalImageType *moving,int level);
												void MoveToNextLevelInThePyramid(InternalImageType *Previous, InternalImageType* NextLevel);

        };//Class Multi Scale Optical Flow

        
} //namespace ITK

#ifndef ITK_MANUAL_INSTANTIATION
#include "itkMultiScaleOpticalFlow.hxx"
#endif

#endif
