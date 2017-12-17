#!/bin/sh
# start up vtk \
exec vtk $0 "$@"

source vtkInclude.tcl
source vtkInt.tcl
source colors.tcl

#-----------------------------------------------------------------------
#	create render windows
#-----------------------------------------------------------------------

vtkRenderMaster master
set ren_win [master MakeRenderWindow]
set ren1 [$ren_win MakeRenderer]
set iren [$ren_win MakeRenderWindowInteractor]

# build up surface

vtkPolyReader surface
  surface SetFilename "surface.dat"
  surface Update

vtkSmoothPolyFilter filter
  filter SetInput [surface GetOutput]

vtkElevationFilter elevation
  elevation SetInput [filter GetOutput]
  elevation SetLowPoint 0 0 0
  elevation SetHighPoint 0 1 0
  elevation Update

vtkCastToConcrete cast
  cast SetInput [elevation GetOutput]
  
vtkLookupTable lookup
  lookup SetHueRange 0.2 1
  lookup SetSaturationRange 0.8 1
  lookup SetValueRange 0.9 0.9
  eval lookup SetTableRange [[elevation GetOutput] GetScalarRange]
  puts [[elevation GetOutput] GetScalarRange]

vtkPolyMapper surf_mapper
  surf_mapper SetInput [cast GetPolyDataOutput]
  surf_mapper SetLookupTable lookup
  eval surf_mapper SetScalarRange [lookup GetTableRange]

vtkActor surf_actor
  surf_actor SetMapper surf_mapper
#  eval [surf_actor GetProperty] SetColor $navy

$ren1 AddActors surf_actor

# build up outline

vtkOutlineFilter outline
  outline SetInput [surface GetOutput]

vtkPolyMapper outline_mapper
  outline_mapper SetInput [outline GetOutput]

vtkActor outline_actor
  outline_actor SetMapper outline_mapper
  eval [outline_actor GetProperty] SetColor $white

$ren1 AddActors outline_actor

# finish initalization

eval $ren1 SetBackground $grey

$iren Initialize
$iren SetUserMethod { wm deiconify .vtkInteract }
$ren_win Render

wm withdraw .
