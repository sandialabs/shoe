/*=========================================================================

  Copyright 2012 Sandia Corporation, Ken Martin, Will Schroeder, Bill Lorensen
  All rights reserved.
  See Copyright.txt or http://www.kitware.com/Copyright.htm for details.

=========================================================================*/
// .NAME vtkShoeBoxAlgorithm - Superclass for algorithms that produce only SHOE higher order meshes as output
// .SECTION Description

// vtkShoeBoxAlgorithm is a convenience class to make writing algorithms
// easier. It is also designed to help transition old algorithms to the new
// pipeline architecture. Ther are some assumptions and defaults made by this
// class you should be aware of. This class defaults such that your filter
// will have one input port and one output port. If that is not the case
// simply change it with SetNumberOfInputPorts etc. See this classes
// constructor for the default. This class also provides a FillInputPortInfo
// method that by default says that all inputs will be ShoeBox. If that
// isn't the case then please override this method in your subclass. This
// class breaks out the downstream requests into seperate functions such as
// ExecuteData and ExecuteInformation.  For new algorithms you should
// implement RequestData( request, inputVec, outputVec) but for older filters
// there is a default implementation that calls the old ExecuteData(output)
// signature, for even older filters that don't implement ExecuteData the
// default implementation calls the even older Execute() signature.

#ifndef __vtkShoeBoxAlgorithm_h
#define __vtkShoeBoxAlgorithm_h

#include "vtkAlgorithm.h"

class vtkDataSet;
class vtkShoeBox;

class VTK_FILTERING_EXPORT vtkShoeBoxAlgorithm : public vtkAlgorithm
{
public:
  static vtkShoeBoxAlgorithm *New();
  vtkTypeRevisionMacro(vtkShoeBoxAlgorithm,vtkAlgorithm);
  void PrintSelf(ostream& os, vtkIndent indent);

  // Description:
  // Get the output data object for a port on this algorithm.
  vtkShoeBox* GetOutput();
  vtkShoeBox* GetOutput(int);
  virtual void SetOutput(vtkDataObject* d);

  // Description:
  // see vtkAlgorithm for details
  virtual int ProcessRequest(vtkInformation*,
                             vtkInformationVector**,
                             vtkInformationVector*);

  // this method is not recommended for use, but lots of old style filters
  // use it
  vtkDataObject *GetInput(int port);
  vtkDataObject *GetInput() { return this->GetInput(0); };
  vtkShoeBox *GetShoeBoxInput(int port);

  // Description:
  // Set an input of this algorithm.
  void SetInput(vtkDataObject *);
  void SetInput(int, vtkDataObject*);

  // Description:
  // Add an input of this algorithm.
  void AddInput(vtkDataObject *);
  void AddInput(int, vtkDataObject*);

protected:
  vtkShoeBoxAlgorithm();
  ~vtkShoeBoxAlgorithm();

  // convenience method
  virtual int RequestInformation(vtkInformation* request,
                                 vtkInformationVector** inputVector,
                                 vtkInformationVector* outputVector);

  // Description:
  // This is called by the superclass.
  // This is the method you should override.
  virtual int RequestData(vtkInformation* request,
                          vtkInformationVector** inputVector,
                          vtkInformationVector* outputVector);

  // Description:
  // This is called by the superclass.
  // This is the method you should override.
  virtual int RequestUpdateExtent(vtkInformation*,
                                  vtkInformationVector**,
                                  vtkInformationVector*);

  // Description:
  // This method is the old style execute method
  virtual void ExecuteData(vtkDataObject *output);
  virtual void Execute();

  // see algorithm for more info
  virtual int FillOutputPortInformation(int port, vtkInformation* info);
  virtual int FillInputPortInformation(int port, vtkInformation* info);

private:
  vtkShoeBoxAlgorithm(const vtkShoeBoxAlgorithm&);  // Not implemented.
  void operator=(const vtkShoeBoxAlgorithm&);  // Not implemented.
};

#endif
