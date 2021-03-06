
LOCAL=/usr/local/lib
CXX = c++
LINK.o = $(CXX) $(LDFLAGS)
CXXFLAGS = -I/usr/local/include/vtk -I../include -Wall -DDEBUG -g -O3
LDFLAGS = -lMesaGL -lVTKCommon -lVTKGraphics -L/usr/X11/lib -lXext \
	-lXt -lX11 -lg++

SOURCES = bspline.cc surfgen.cc real2spl.cc approximation.cc \
	approximation1.cc approximation2.cc approximation3.cc \
	approximation4.cc householder.cc functional.cc \
	polygon.cc wave.cc set.cc exclusion.cc spl2vtk.cc spl2real.cc \
	cmpspl.cc

all:	spl2real real2spl spl2vtk cmpspl

sample: exclusion.o set.o polygon.o
	$(CXX) -o $@ $^ $(LDFLAGS)

wave: wave.o 
	$(CXX) -o $@ $^ $(LDFLAGS)
	ln -f wave ../

real2spl: real2spl.o approximation.o approximation1.o approximation2.o \
	approximation3.o approximation4.o functional.o householder.o \
	bspline.o exclusion.o set.o polygon.o
	$(CXX) -o $@ $^ $(LDFLAGS)
	ln -f real2spl ../

spl2real: spl2real.o surfgen.o bspline.o exclusion.o set.o polygon.o
	$(CXX) -o $@ $^ $(LDFLAGS)
	ln -f spl2real ../

spl2vtk: spl2vtk.o surfgen.o bspline.o exclusion.o set.o polygon.o
	$(CXX) -o $@ $^ $(LDFLAGS)
	ln -f spl2vtk ../

cmpspl: cmpspl.o bspline.o
	$(CXX) -o $@ $^ $(LDFLAGS)
	ln -f cmpspl ../

backup:
	cd ../../; find bspline \( -name "*.o" -or -name "*~" \) \
	 -print >/tmp/$$$$; \
	tar cvzf $(HOME)/backup/bspline.tar.gz -X /tmp/$$$$ bspline ;\
	rm /tmp/$$$$

depend:
	makedepend -Y $(CXXFLAGS) $(SOURCES)

clean:
	rm *.o core *~ spl2real real2spl spl2vtk

# DO NOT DELETE THIS LINE -- make depend depends on it.

bspline.o: ../include/bspline.hh ../include/stdtypes.h ../include/array.hh
bspline.o: ../include/debug.hh ../include/array.cc ../include/funcs.hh
bspline.o: ../include/list.hh ../include/list.cc ../include/vector.hh
bspline.o: ../include/vector.cc ../include/angles.hh ../include/map.hh
bspline.o: ../include/btree.hh ../include/pair.hh ../include/btree.cc
bspline.o: ../include/map.cc ../include/matrix.hh ../include/matrix.cc
surfgen.o: ../include/surfgen.hh ../include/bspline.hh ../include/stdtypes.h
surfgen.o: ../include/array.hh ../include/debug.hh ../include/array.cc
surfgen.o: ../include/funcs.hh ../include/list.hh ../include/list.cc
surfgen.o: ../include/vector.hh ../include/vector.cc ../include/angles.hh
surfgen.o: ../include/map.hh ../include/btree.hh ../include/pair.hh
surfgen.o: ../include/btree.cc ../include/map.cc ../include/matrix.hh
surfgen.o: ../include/matrix.cc /usr/local/include/vtk/vtkFloatPoints.h
surfgen.o: /usr/local/include/vtk/vtkPoints.h
surfgen.o: /usr/local/include/vtk/vtkRefCount.h
surfgen.o: /usr/local/include/vtk/vtkObject.h
surfgen.o: /usr/local/include/vtk/vtkTimeStamp.h
surfgen.o: /usr/local/include/vtk/vtkSetGet.h
surfgen.o: /usr/local/include/vtk/vtkIndent.h
surfgen.o: /usr/local/include/vtk/vtkIdList.h
surfgen.o: /usr/local/include/vtk/vtkIntArray.h
surfgen.o: /usr/local/include/vtk/vtkFloatArray.h
surfgen.o: /usr/local/include/vtk/vtkFloatScalars.h
surfgen.o: /usr/local/include/vtk/vtkScalars.h
surfgen.o: /usr/local/include/vtk/vtkCellArray.h
surfgen.o: /usr/local/include/vtk/vtkCell.h
surfgen.o: /usr/local/include/vtk/vtkCellType.h
surfgen.o: /usr/local/include/vtk/vtkPolyData.h
surfgen.o: /usr/local/include/vtk/vtkPointSet.h
surfgen.o: /usr/local/include/vtk/vtkDataSet.h
surfgen.o: /usr/local/include/vtk/vtkPointData.h
surfgen.o: /usr/local/include/vtk/vtkVectors.h
surfgen.o: /usr/local/include/vtk/vtkFloatVectors.h
surfgen.o: /usr/local/include/vtk/vtkNormals.h
surfgen.o: /usr/local/include/vtk/vtkFloatNormals.h
surfgen.o: /usr/local/include/vtk/vtkTCoords.h
surfgen.o: /usr/local/include/vtk/vtkFloatTCoords.h
surfgen.o: /usr/local/include/vtk/vtkTensors.h
surfgen.o: /usr/local/include/vtk/vtkTensor.h
surfgen.o: /usr/local/include/vtk/vtkFloatTensors.h
surfgen.o: /usr/local/include/vtk/vtkUserDefined.h
surfgen.o: /usr/local/include/vtk/vtkVoidArray.h
surfgen.o: /usr/local/include/vtk/vtkPointLocator.h
surfgen.o: /usr/local/include/vtk/vtkLocator.h
surfgen.o: /usr/local/include/vtk/vtkCellList.h
surfgen.o: /usr/local/include/vtk/vtkLinkList.h
surfgen.o: /usr/local/include/vtk/vtkStructuredGrid.h
surfgen.o: /usr/local/include/vtk/vtkStructuredData.h
surfgen.o: /usr/local/include/vtk/vtkBitScalars.h
surfgen.o: /usr/local/include/vtk/vtkBitArray.h
surfgen.o: /usr/local/include/graphics/vtkPolyWriter.h
surfgen.o: /usr/local/include/graphics/vtkDataWriter.h
surfgen.o: /usr/local/include/graphics/vtkWriter.h
surfgen.o: /usr/local/include/graphics/vtkStructuredGridWriter.h
surfgen.o: ../include/exclusion.hh ../include/set.hh ../include/polygon.hh
real2spl.o: ../include/approximation.hh ../include/bspline.hh
real2spl.o: ../include/stdtypes.h ../include/array.hh ../include/debug.hh
real2spl.o: ../include/array.cc ../include/funcs.hh ../include/list.hh
real2spl.o: ../include/list.cc ../include/vector.hh ../include/vector.cc
real2spl.o: ../include/angles.hh ../include/map.hh ../include/btree.hh
real2spl.o: ../include/pair.hh ../include/btree.cc ../include/map.cc
real2spl.o: ../include/matrix.hh ../include/matrix.cc ../include/exclusion.hh
real2spl.o: ../include/set.hh ../include/polygon.hh ../include/functional.hh
approximation.o: ../include/approximation.hh ../include/bspline.hh
approximation.o: ../include/stdtypes.h ../include/array.hh
approximation.o: ../include/debug.hh ../include/array.cc ../include/funcs.hh
approximation.o: ../include/list.hh ../include/list.cc ../include/vector.hh
approximation.o: ../include/vector.cc ../include/angles.hh ../include/map.hh
approximation.o: ../include/btree.hh ../include/pair.hh ../include/btree.cc
approximation.o: ../include/map.cc ../include/matrix.hh ../include/matrix.cc
approximation.o: ../include/exclusion.hh ../include/set.hh
approximation.o: ../include/polygon.hh ../include/functional.hh
approximation.o: ../include/householder.hh /usr/local/include/vtk/vtkMath.h
approximation1.o: ../include/approximation.hh ../include/bspline.hh
approximation1.o: ../include/stdtypes.h ../include/array.hh
approximation1.o: ../include/debug.hh ../include/array.cc ../include/funcs.hh
approximation1.o: ../include/list.hh ../include/list.cc ../include/vector.hh
approximation1.o: ../include/vector.cc ../include/angles.hh ../include/map.hh
approximation1.o: ../include/btree.hh ../include/pair.hh ../include/btree.cc
approximation1.o: ../include/map.cc ../include/matrix.hh ../include/matrix.cc
approximation1.o: ../include/exclusion.hh ../include/set.hh
approximation1.o: ../include/polygon.hh ../include/functional.hh
approximation1.o: ../include/householder.hh /usr/local/include/vtk/vtkMath.h
approximation2.o: ../include/approximation.hh ../include/bspline.hh
approximation2.o: ../include/stdtypes.h ../include/array.hh
approximation2.o: ../include/debug.hh ../include/array.cc ../include/funcs.hh
approximation2.o: ../include/list.hh ../include/list.cc ../include/vector.hh
approximation2.o: ../include/vector.cc ../include/angles.hh ../include/map.hh
approximation2.o: ../include/btree.hh ../include/pair.hh ../include/btree.cc
approximation2.o: ../include/map.cc ../include/matrix.hh ../include/matrix.cc
approximation2.o: ../include/exclusion.hh ../include/set.hh
approximation2.o: ../include/polygon.hh ../include/functional.hh
approximation2.o: ../include/householder.hh /usr/local/include/vtk/vtkMath.h
approximation3.o: ../include/approximation.hh ../include/bspline.hh
approximation3.o: ../include/stdtypes.h ../include/array.hh
approximation3.o: ../include/debug.hh ../include/array.cc ../include/funcs.hh
approximation3.o: ../include/list.hh ../include/list.cc ../include/vector.hh
approximation3.o: ../include/vector.cc ../include/angles.hh ../include/map.hh
approximation3.o: ../include/btree.hh ../include/pair.hh ../include/btree.cc
approximation3.o: ../include/map.cc ../include/matrix.hh ../include/matrix.cc
approximation3.o: ../include/exclusion.hh ../include/set.hh
approximation3.o: ../include/polygon.hh ../include/functional.hh
approximation3.o: ../include/householder.hh /usr/local/include/vtk/vtkMath.h
approximation4.o: ../include/approximation.hh ../include/bspline.hh
approximation4.o: ../include/stdtypes.h ../include/array.hh
approximation4.o: ../include/debug.hh ../include/array.cc ../include/funcs.hh
approximation4.o: ../include/list.hh ../include/list.cc ../include/vector.hh
approximation4.o: ../include/vector.cc ../include/angles.hh ../include/map.hh
approximation4.o: ../include/btree.hh ../include/pair.hh ../include/btree.cc
approximation4.o: ../include/map.cc ../include/matrix.hh ../include/matrix.cc
approximation4.o: ../include/exclusion.hh ../include/set.hh
approximation4.o: ../include/polygon.hh ../include/functional.hh
approximation4.o: ../include/householder.hh /usr/local/include/vtk/vtkMath.h
householder.o: ../include/householder.hh ../include/matrix.hh
householder.o: ../include/debug.hh ../include/stdtypes.h ../include/matrix.cc
functional.o: ../include/functional.hh ../include/matrix.hh
functional.o: ../include/debug.hh ../include/stdtypes.h ../include/matrix.cc
functional.o: ../include/bspline.hh ../include/array.hh ../include/array.cc
functional.o: ../include/funcs.hh ../include/list.hh ../include/list.cc
functional.o: ../include/vector.hh ../include/vector.cc ../include/angles.hh
functional.o: ../include/map.hh ../include/btree.hh ../include/pair.hh
functional.o: ../include/btree.cc ../include/map.cc
polygon.o: ../include/polygon.hh ../include/debug.hh ../include/stdtypes.h
polygon.o: ../include/array.hh ../include/array.cc ../include/funcs.hh
polygon.o: /usr/local/include/vtk/vtkMath.h
wave.o: ../include/matrix.hh ../include/debug.hh ../include/stdtypes.h
wave.o: ../include/matrix.cc ../include/bspline.hh ../include/array.hh
wave.o: ../include/array.cc ../include/funcs.hh ../include/list.hh
wave.o: ../include/list.cc ../include/vector.hh ../include/vector.cc
wave.o: ../include/angles.hh ../include/map.hh ../include/btree.hh
wave.o: ../include/pair.hh ../include/btree.cc ../include/map.cc
set.o: ../include/set.hh ../include/debug.hh ../include/stdtypes.h
set.o: ../include/funcs.hh
exclusion.o: ../include/exclusion.hh ../include/debug.hh
exclusion.o: ../include/stdtypes.h ../include/list.hh ../include/list.cc
exclusion.o: ../include/funcs.hh ../include/set.hh ../include/polygon.hh
exclusion.o: ../include/array.hh ../include/array.cc
spl2vtk.o: ../include/surfgen.hh ../include/bspline.hh ../include/stdtypes.h
spl2vtk.o: ../include/array.hh ../include/debug.hh ../include/array.cc
spl2vtk.o: ../include/funcs.hh ../include/list.hh ../include/list.cc
spl2vtk.o: ../include/vector.hh ../include/vector.cc ../include/angles.hh
spl2vtk.o: ../include/map.hh ../include/btree.hh ../include/pair.hh
spl2vtk.o: ../include/btree.cc ../include/map.cc ../include/matrix.hh
spl2vtk.o: ../include/matrix.cc ../include/exclusion.hh ../include/set.hh
spl2vtk.o: ../include/polygon.hh
spl2real.o: ../include/surfgen.hh ../include/bspline.hh ../include/stdtypes.h
spl2real.o: ../include/array.hh ../include/debug.hh ../include/array.cc
spl2real.o: ../include/funcs.hh ../include/list.hh ../include/list.cc
spl2real.o: ../include/vector.hh ../include/vector.cc ../include/angles.hh
spl2real.o: ../include/map.hh ../include/btree.hh ../include/pair.hh
spl2real.o: ../include/btree.cc ../include/map.cc ../include/matrix.hh
spl2real.o: ../include/matrix.cc
cmpspl.o: ../include/bspline.hh ../include/stdtypes.h ../include/array.hh
cmpspl.o: ../include/debug.hh ../include/array.cc ../include/funcs.hh
cmpspl.o: ../include/list.hh ../include/list.cc ../include/vector.hh
cmpspl.o: ../include/vector.cc ../include/angles.hh ../include/map.hh
cmpspl.o: ../include/btree.hh ../include/pair.hh ../include/btree.cc
cmpspl.o: ../include/map.cc ../include/matrix.hh ../include/matrix.cc
