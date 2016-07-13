TEMPLATE = subdirs

SUBDIRS += \
    Denoising \ 
    OpenMesh     

Denoising.depends = OpenMesh
