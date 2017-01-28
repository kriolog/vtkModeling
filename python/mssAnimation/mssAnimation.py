import sys
from vtk import *
import Tkinter
from vtk.tk.vtkTkRenderWindowInteractor import \
     vtkTkRenderWindowInteractor

if len(sys.argv) < 3:
    raise Exception, "Usage: python mssAnimation.py min max"

mini = int(sys.argv[1])
maxi = int(sys.argv[2])
current = mini

reader=vtkDataSetReader()
reader.SetFileName(str(mini))
reader.Modified()
reader.Update()
mapper = vtkDataSetMapper()
mapper.SetInput(reader.GetOutput())
actor = vtkActor()
actor.SetMapper(mapper)
ren = vtkRenderer()
ren.AddActor(actor)
ren.SetBackground(0.0,0.0,0.5)
renWin = vtkRenderWindow()
renWin.AddRenderer(ren)

def PreviousImage():
    global current
    if current > mini:
        reader.SetFileName(str(current))
        current = current - 1
        reader.Modified()
        reader.Update()
        renWin. Render()

def NextImage():
    global current
    if maxi > current:
        reader.SetFileName(str(current))
        current = current + 1
        reader.Modified()
        reader.Update()
        renWin. Render()

def Display():
    global current
    for i in range(mini,maxi*2):
        if current == maxi:
            current = 0
        else:
            current = current + 1
        reader.SetFileName(str(current))
        reader.Modified()
        reader.Update()
        renWin. Render()       

# Now actually create the GUI
root = Tkinter.Tk()
root.withdraw()
top = Tkinter.Toplevel(root)

# Define a quit method that exits cleanly.
def quit(obj=root):
    obj.quit()

display_frame = Tkinter.Frame(top)
display_frame.pack(side="top", anchor="n", fill="both", expand="false")

# Buttons
ctrl_buttons = Tkinter.Frame(top)
ctrl_buttons.pack(side="top", anchor="n", fill="both", expand="false")
quit_button = Tkinter.Button(ctrl_buttons, text="Quit", command=quit)
next_button = Tkinter.Button(ctrl_buttons, text="Next",
                                command=NextImage)
previous_button = Tkinter.Button(ctrl_buttons, text="Previous",
                                command=PreviousImage)
display_button = Tkinter.Button(ctrl_buttons, text="Movie",
                                command=Display)

for i in (quit_button, previous_button, next_button, display_button):
    i.pack(side="left", expand="true", fill="both")


# Create the render widget
renderer_frame = Tkinter.Frame(display_frame)
renderer_frame.pack(padx=3, pady=3,side="left", anchor="n",
                    fill="both", expand="false")

render_widget = vtkTkRenderWindowInteractor(renderer_frame,
                                            rw=renWin, width=800,
                                            height=800)
render_widget.GetRenderWindow().GetInteractor().SetInteractorStyle(vtkInteractorStyleTrackballCamera())
for i in (render_widget, display_frame):
    i.pack(side="top", anchor="n",fill="both", expand="false")


# Done with the GUI. 
###
iact = render_widget.GetRenderWindow().GetInteractor()

# Create an initial interesting view
cam1 = ren.GetActiveCamera()
cam1.Elevation(110)
cam1.SetViewUp(0, 0, -1)
cam1.Azimuth(45)
ren.ResetCameraClippingRange()

# Render it
render_widget.Render()

iact.Initialize()
renWin.Render()
iact.Start()

# Start Tkinter event loop
root.mainloop()
