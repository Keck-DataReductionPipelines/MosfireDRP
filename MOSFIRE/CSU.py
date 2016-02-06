

''' 
MOSFIRE CSU Utlity Code
    Includes physical parameters of CSU.

Created March 2, 2011 by npk

numslits, numbars give information about the number of slits and bars in the CSU.

Note, for consistency with the hardware bars and slits are indexed from 1.

tempscale is the thermal scaling factor to shrink room temperature linear 
dimensions to 120 K. The number 0.99646 is from R. Weber's spreadsheet 
"MOSFIRE Thermal Dimension Scaling Factors.xls", referenced from 
"Thermophysical properties of matter, Vol 12-13, Thermal Expansion".

demagnification (7.24254), and center_pix (1042.99 pix, 1035.88 pix) is 
measured by ccs in January and Feb 2011 using pinhole mask data taken 
during the eighth cooldown. These are described in 
"README.focal_plane_mapping.txt".

bar pitch is 5.8 mm which is related to pixels using the demagnification 
and temperature scale.

'''
import Detector, IO
import numpy as np
import unittest
import os
#from pyraf import iraf

import pdb
from MosfireDrpLog import debug, info, warning, error

from astropy.modeling import models
import pickle

class MismatchError(Exception):
    '''The code expected a CSU with 46 slits, but found something else.'''

    def __init__(self, value):
        self.parameter = value

numslits = 46
numbars = numslits * 2

tempscale = 0.99646 

def mm_to_pix(mm):
    return mm/demagnification/Detector.pixelsize

mm = 1
demagnification = 7.24254
center_pix = (1042.986, 1035.879)
barpitch_mm = (5.8 * mm * tempscale)
barpitch_pix = mm_to_pix(5.8 * mm * tempscale)

def in_field(px, py):
    '''Determines if the pixel coordinate (x,y) is within the circular 1150 pix FOV'''

    x = px-center_pix[0]
    y = py-center_pix[1]
    dist = np.sqrt(x*x + y*y)
    if dist < 1150.: return True
    return False

def csu_pix_to_mm_poly(x_pix, y_pix):
    #(x_kfp, y_kfp) = mosfire_geoxytran(x_pix, y_pix, direction="backward")
    (x_kfp, y_kfp) = python_geoxytran(x_pix, y_pix, direction="backward")
    centerx = 137.400

    x_mm = centerx - x_kfp
    y_mm = y_kfp

    return (x_mm, y_mm)

def csu_mm_to_pix_poly(x_mm, slitno):
    '''Uses ccs fits in ../platescale directory'''
    # _kfp is keck focal plane
    centerx = 137.400
    x_kfp = (centerx - x_mm) 
    y_kfp = 5.8 * (numslits/2. - slitno + 0.35)  * tempscale


    #return mosfire_geoxytran(x_kfp, y_kfp)
    return python_geoxytran(x_kfp, y_kfp)


# def csu_mm_to_pix(x_mm, slitno, Please_Use=False):
#     '''Convert a slit's position into a pixel value. This is a linear approximation to a sixth order polynomial fit by ccs.
#     Positions are index from 1: 1 .. 2048
#     '''

#     if Please_Use==False:
#         error("Use csu_mm_to_pix_poly (a polynomial fit) rather than csu_mm_to_pix (a linear fit)")
#         raise Exception("Use csu_mm_to_pix_poly (a polynomial fit) rather than csu_mm_to_pix (a linear fit)")
#         return

#     # _kfp is keck focal plane
#     centerx = 137.400
#     x_kfp = (centerx - x_mm) * tempscale
#     y_kfp = 5.8*mm * (numslits/2. - slitno + 0.35) * tempscale

#     path = os.path.join(os.path.dirname(__file__), "data", 
#             "linear_pix2mm_120k.db") 
#     #
#     return  mosfire_geoxytran(x_kfp, y_kfp,
#             database=path, transform="linear_pix2mm_120k")


def mosfire_geoxytran(x_kfp, y_kfp, transform="final.pix2mm.4.972.120k",
        database="/platescale/10March2011.4.972.db", direction="forward"):
    '''Conveninece wrapper around IRAF geoxytran'''
    iraf.images()


    path = os.path.join(os.path.dirname(__file__), "data", 
            "10March2011.4.972.db")
    database = path
    pars = iraf.geoxytran.getParList()
    iraf.geoxytran.unlearn()
    t = iraf.geoxytran("STDIN", "STDOUT", Stdin=["%f %f" % (x_kfp, y_kfp)], Stdout=1,
            database=database,
            transform=transform,
            direction=direction)

    iraf.geoxytran.setParList(pars) 
    (x,y) = np.array(t[0].split(), dtype=np.float64)

    return (x,y)

def python_geoxytran(x_kfp, y_kfp, direction="forward"):
    if direction=="backward":
        px_backward = models.Legendre2D(6, 6, c0_0=-1.359159118077272, c1_0=0.12849509682026816, c2_0=0.00017587310282272408, c3_0=-8.214009406649863e-06,
                                   c4_0=2.0206624190921399e-07, c5_0=-2.225331028379213e-09, c6_0=1.8201097072390407e-14, c0_1=-0.0008865780983085235,
                                   c1_1=0.0002837258293901996, c2_1=-1.3953479926814954e-05, c3_1=2.6865725414225316e-07, c4_1=3.7333388292351965e-10,
                                   c5_1=8.662694037494459e-13, c6_1=-7.6802320344598e-15, c0_2=9.616027720780746e-05, c1_2=-1.5463196269995818e-05,
                                   c2_2=4.615248418093103e-07, c3_2=-9.938940430240368e-09, c4_2=7.538338883385691e-12, c5_2=-1.8016883087953452e-13,
                                   c6_2=1.6284543459178821e-15, c0_3=-2.7112157097925163e-06, c1_3=3.63691974893539e-07, c2_3=6.326334170834919e-10,
                                   c3_3=1.2045620539279474e-11, c4_3=-6.281301326529564e-13, c5_3=1.5395969945758583e-14, c6_3=-1.4203191615580046e-16,
                                   c0_4=2.487234831550635e-08, c1_4=-5.3202681529753e-09, c2_4=3.813876920246599e-12, c3_4=-4.578771786695712e-13,
                                   c4_4=2.4833675429790513e-14, c5_4=-6.278532214053127e-16, c6_4=5.932362209122972e-18, c0_5=2.6533817724685113e-10,
                                   c1_5=6.362774492493808e-14, c2_5=-5.695287662674329e-14, c3_5=7.648943667217284e-15, c4_5=-4.4244874441233506e-16,
                                   c5_5=1.1718033882619874e-17, c6_5=-1.1450561454795142e-19, c0_6=-5.252495563272626e-15, c1_6=6.498737275590606e-16,
                                   c2_6=2.2221508682832634e-16, c3_6=-4.096197448486931e-17, c4_6=2.7086424901520096e-18, c5_6=-7.787892566015997e-20,
                                   c6_6=8.028451974197805e-22)
        py_backward = models.Legendre2D(6, 6, c0_0=-1.3408760539296245, c1_0=-0.0014681933080717899, c2_0=6.252434078059442e-05, c3_0=-1.7794023960848901e-06,
                                   c4_0=2.0505693079301286e-08, c5_0=1.3303121908968087e-10, c6_0=4.036925907590215e-14, c0_1=0.1287659978047137,
                                   c1_1=0.0002187658143857909, c2_1=-1.1414122040749694e-05, c3_1=2.514881941931133e-07, c4_1=-4.014646650126551e-09,
                                   c5_1=4.6361664655461665e-12, c6_1=-4.2907954493018394e-14, c0_2=0.00017210509816287917, c1_2=-1.1517572721650909e-05,
                                    c2_2=4.070900780580943e-07, c3_2=-2.924032092730881e-10, c4_2=4.3651272195759074e-11, c5_2=-1.0942185864222553e-12,
                                    c6_2=1.0124374619198603e-14, c0_3=-8.927312822514692e-06, c1_3=2.178919731355544e-07, c2_3=-9.309247977587529e-09,
                                    c3_3=7.587241655284752e-11, c4_3=-4.20296301764774e-12, c5_3=1.0534116340750084e-13, c6_3=-9.745842383678927e-16,
                                    c0_4=2.3500427475161216e-07, c1_4=5.797618861500032e-10, c2_4=2.787756737792332e-11, c3_4=-3.471711550473785e-12,
                                     c4_4=1.9236825440442053e-13, c5_4=-4.82226129561789e-15, c6_4=4.4622762166743036e-17, c0_5=-2.631606790683655e-09,
                                     c1_5=1.8088697785601813e-12, c2_5=-6.029913349859671e-13, c3_5=7.51954755289406e-14, c4_5=-4.1690939247203245e-15,
                                     c5_5=1.0455627246423308e-16, c6_5=-9.679289662228309e-19, c0_6=1.1235020215574227e-14, c1_6=-1.481115941278654e-14,
                                     c2_6=4.964545148330356e-15, c3_6=-6.201605688722248e-16, c4_6=3.441248923238193e-17, c5_6=-8.635683356739731e-19, c6_6=7.999236760366155e-21)
        x = px_backward(x_kfp/100.0,y_kfp/100.0)*100.0
        y = py_backward(x_kfp/100.0,y_kfp/100.0)*100.0

    
    if direction=="forward":
        px_forward = models.Legendre2D(6, 6, c0_0=10.429995608040636, c1_0=7.669072866696911, c2_0=-0.0005171642200091058, c3_0=0.0010640118370285666,
                                    c4_0=6.591355825164487e-05, c5_0=0.0004579785406086718, c6_0=-3.381749890396372e-07, c0_1=-0.029899970699651657,
                                    c1_1=-1.0370107736602057e-05, c2_1=0.00012355149488906862, c3_1=0.0001681000870132385, c4_1=-8.915078035548195e-05,
                                    c5_1=1.0828480948981245e-06, c6_1=-3.0604028638458465e-07, c0_2=-0.0005617773576709009, c1_2=0.0021497066157591966,
                                    c2_2=0.0003159972245946561, c3_2=0.001999515078707485, c4_2=9.953809608627005e-07, c5_2=8.667245967324062e-06,
                                    c6_2=2.8043195865254592e-08, c0_3=0.0001440482443085916, c1_3=-3.4260998883389794e-05, c2_3=-0.00016054621466681272,
                                    c3_3=9.517312759587115e-07, c4_3=-4.1446577769773705e-07, c5_3=3.6262929000660604e-07, c6_3=-1.4504338440561597e-07,
                                    c0_4=0.00015388407922332725, c1_4=0.0010841802946648064, c2_4=5.000639720902258e-07, c3_4=5.927690431156899e-06,
                                    c4_4=5.991902315389979e-07, c5_4=1.7610373474546858e-06, c6_4=3.6698527469006115e-09, c0_5=-5.883225927240384e-05,
                                    c1_5=3.472735736376934e-07, c2_5=-1.0563775236811306e-06, c3_5=2.2897505242989876e-07, c4_5=-2.389407023502443e-07,
                                    c5_5=-2.7244280935829703e-09, c6_5=-2.0771844124138064e-08, c0_6=2.8762976671867224e-07, c1_6=3.2864844277261225e-06,
                                    c2_6=3.881850299314179e-07, c3_6=1.3047967293000456e-06, c4_6=1.0513711949538813e-08, c5_6=8.671289818897808e-09, c6_6=-1.760160785036003e-09)
        py_forward = models.Legendre2D(6, 6, c0_0=10.35865352089189, c1_0=0.029868752847050786, c2_0=-5.181441268292675e-05, c3_0=9.57356706812987e-05,
                                    c4_0=0.00017241278829382978, c5_0=-2.4675256900133155e-05, c6_0=3.95235509529886e-07, c0_1=7.66897227133641,
                                    c1_1=0.00034364965616906756, c2_1=0.0027356854767171535, c3_1=0.00042743941177974913, c4_1=0.0007837226536120123,
                                    c5_1=9.874568475762985e-07, c6_1=2.3707441689970604e-06, c0_2=0.0003489725952788292, c1_2=0.0002812999003417556,
                                    c2_2=0.0006371193460650293, c3_2=-8.6155293329177e-05, c4_2=1.6936419268906433e-06, c5_2=-3.490309288135124e-07,
                                    c6_2=4.3797605194396266e-07, c0_3=0.0007884544919103953, c1_3=0.00021149527310720538, c2_3=0.0017752319143687079,
                                    c3_3=1.4502993553656453e-06, c4_3=5.374215134624355e-06, c5_3=4.6355913497424066e-07, c6_3=9.418886502384274e-07,
                                    c0_4=9.01336941821207e-05, c1_4=-0.0001230909091922516, c2_4=2.2577487548726405e-06, c3_4=8.480105722621575e-08,
                                     c4_4=6.821359101651797e-07, c5_4=-2.7924348331191066e-07, c6_4=1.195944732829135e-08, c0_5=0.0005414644597315499,
                                     c1_5=1.1375883331712563e-06, c2_5=8.954534272097303e-06, c3_5=5.176871464493416e-07, c4_5=1.5475712624698004e-06,
                                     c5_5=1.0885690389676392e-08, c6_5=1.1381062005799344e-08, c0_6=1.4659708681089706e-07, c1_6=-6.338488221563267e-07,
                                     c2_6=2.1673413055835894e-07, c3_6=-2.8336747286630153e-07, c4_6=1.6623075723212755e-08, c5_6=1.7700444635232546e-08, c6_6=-1.7039889450896992e-09)
        x = px_forward(x_kfp/100.0,y_kfp/100.0)*100.0
        y = py_forward(x_kfp/100.0,y_kfp/100.0)*100.0



    return (x,y)


def mosfire_geoxytrans(x_kfp, y_kfp, transform="final.pix2mm.4.972.120k",
        database="ale/10March2011.4.972.db", direction="forward"):
    '''Conveninece wrapper around IRAF geoxytran'''
    iraf.images()
    path = os.path.join(os.path.dirname(__file__), "data", 
            "10March2011.4.972.db")
    database = path
    pars = iraf.geoxytran.getParList()
    iraf.geoxytran.unlearn()
    ins = []
    for i in range(len(x_kfp)):
        ins.append("%f %f" % (x_kfp[i], y_kfp[i]))
    results = iraf.geoxytran("STDIN", "STDOUT", Stdin=ins, Stdout=1,
            database=database,
            transform=transform,
            direction=direction)

    iraf.geoxytran.setParList(pars) 

    poss = []
    for result in results:
        poss.append(np.array(result.split(), dtype=np.float64))

    return np.array(poss)

    
def bar_to_slit(x):
    '''Convert a bar #(1-92) to a slit(1-46) number'''
    if (x < 1) or (x > numbars):
        error("Not indexing CSU properly")
        raise MismatchError("Not indexing CSU properly")
    return int(x+1)/2

def to_ds9_region(poss, dash=1, color="green", label=True):
    s = []
    d = np.radians(4.2)
    dx = barpitch_pix/2. * np.sin(d)
    dy = barpitch_pix/2. * np.cos(d)

    for i in range(1,numbars+1):
        pos = poss[i-1]
        if not np.isfinite(pos[0]): continue
        if not np.isfinite(pos[1]): continue

        ln = [pos[0]+dx,pos[1]-dy,pos[0]-dx,pos[1]+dy]
        if label:
            s.append("line(%6.3f, %6.3f, %6.3f, %6.3f) # line=0 0 color=%s text={b%2.0i} dash=%1i fixed=1 edit=0 move=0 rotate=0 \n" % (ln[0], ln[1], ln[2], ln[3], color, i, dash))
        else:
            s.append("line(%6.3f, %6.3f, %6.3f, %6.3f) # line=0 0 color=%s dash=%1i fixed=1 edit=0 move=0 rotate=0 \n" % (ln[0], ln[1], ln[2], ln[3], color, dash))

    return s


class Barset:
    '''Barset provides convenience functions around a CSU slitmask'''

    pos = [] 
    pos_pix = []

    header = None
    # Science slit list, mechanical slit list, & alignment slit list.
    ssl = None
    msl = None
    asl = None
    targs = None

    long_slit = False
    long2pos_slit = False

    scislit_to_slit = []
    alignment_slits = []


    def __init__(self):
        pass

    def set_header(self, header, ssl=None, msl=None, asl=None, targs=None):
        '''Passed "header" a FITS header dictionary and converts to a Barset'''
        self.pos = np.array(IO.parse_header_for_bars(header))
        self.set_pos_pix()

        self.ssl = ssl
        self.msl = msl
        self.asl = asl
        self.targs = targs

        def is_alignment_slit(slit):
            return (np.float(slit["Target_Priority"]) < 0)

        # If len(ssl) == 0 then the header is for a long slit
        if (header['MASKNAME'] == 'long2pos'):
            info("long2pos mode in CSU slit determination")
            self.long2pos_slit = True

        if (len(ssl) == 0):
        
            self.long_slit = True

            start = np.int(msl[0]["Slit_Number"])
            stop = np.int(msl[-1]["Slit_Number"])


            for mech_slit in msl:
                mech_slit["Target_in_Slit"] = "long"

            self.ssl = np.array([("1", "??", "??", "??", "??", "??", "??", msl[0]['Slit_width'],
                (stop-start+1)*7.6, "0", "long", "0")],
                dtype= [ ('Slit_Number', '|S2'), 
                ('Slit_RA_Hours', '|S2'), ('Slit_RA_Minutes', '|S2'), ('Slit_RA_Seconds', '|S5'),
                ('Slit_Dec_Degrees', '|S3'), ('Slit_Dec_Minutes', '|S2'), ('Slit_Dec_Seconds', '|S5'), 
                ('Slit_width', '|S5'), ('Slit_length', '|S5'), ('Target_to_center_of_slit_distance', '|S5'), 
                ('Target_Name', '|S80'), ('Target_Priority', '|S1')])
            self.scislit_to_slit = [ np.arange(start,stop) ]
            ssl = None

        # Create a map between scislit number and mechanical slit
        # recall that slits count from 1
        if ssl is not None:
            prev = self.msl[0]["Target_in_Slit"]

            v = []

            for science_slit in ssl:
                targ = science_slit["Target_Name"]
                v.append([int(x) for x in self.msl.field("Slit_Number")[np.where(self.msl.field("Target_in_Slit").rstrip() == targ)[0]]])
            self.scislit_to_slit = v

            if (len(self.scislit_to_slit) != len(ssl)) and not (self.long_slit
                    and len(self.scislit_to_slit) == 1):
                error("SSL should match targets in slit")
                raise Exception("SSL should match targets in slit")


    def is_alignment_slitno(self, slitno):
        return (slitno in self.alignment_slits)

    def csu_slit_center(self, slitno):
        '''Returns the mechanical (middle) position of a csu slit in mm'''

        if (slitno < 1) or (slitno > 46):
            error("The requested slit number (%i) does not exist" % 
                    slitno)
            raise Exception("The requested slit number (%i) does not exist" % 
                    slitno)

        os = self.pos[slitno*2 - 2]
        es = self.pos[slitno*2 - 1]

        return (os+es)/2.
        
    def scislit_to_csuslit(self, scislit):
        '''Convert a science slit number to a mechanical slit list'''
        if (scislit < 1) or (scislit > len(self.ssl)+1):
            error("The requested slit number (%i) does not exist" %
                    scislit)
            raise Exception("The requested slit number (%i) does not exist" %
                    scislit)

        return self.scislit_to_slit[scislit-1]
    
    def csu_slit_to_pixel(self, slit):
        '''Convert a CSU slit number to spatial pixel'''
        y0 = 2013

        if (slit < 1) or (slit > 46):
            error("The requested slit number (%i) does not exist" %
                    slit)
            raise Exception("The requested slit number (%i) does not exist" %
                    slit)

        pixel = np.int(y0 - (slit -1) * 44.22)
        return pixel

    def science_slit_to_pixel(self, scislit):
        '''Convert a science slit number to spatial pixel'''

        if (scislit < 1) or (scislit > len(self.ssl)):
            error("The requested science slit number %i does not exist" \
                    % scislit)
            raise Exception("The requested science slit number %i does not exist" \
                    % scislit)

        slits = self.scislit_to_csuslit(scislit)
        debug(str(slits))
        return self.csu_slit_to_pixel(np.median(slits))

    def set_pos_pix(self):
        # _kfp is keck focal plane
        centerx = 137.400
        x_kfp = (centerx - self.pos) 
        slitno = np.ceil(np.arange(1, numbars+1)/2.)
        y_kfp = 5.8 * (numslits/2. - slitno + 0.35)  * tempscale

        #self.pos_pix = mosfire_geoxytrans(x_kfp, y_kfp)
        #self.pos_pix = python_geoxytran(x_kfp, y_kfp)
        tmp_results = python_geoxytran(x_kfp, y_kfp)
        self.pos_pix = np.asarray([list(p) for p in zip(tmp_results[0],tmp_results[1])])

    def to_ds9_region(self):
        poss = []

        for i in range(1,numbars+1):
            poss.append(self.get_bar_pix(i))

        return to_ds9_region(poss)

    def get_bar_pix(self, bar):
        '''Return the pixel position of bar(1-92)'''
        return self.pos_pix[bar-1]


class TestCSUFunctions(unittest.TestCase):

    def setUp(self):
        pass

    def test_bar_to_slit(self):
        sa = self.assertTrue
        sa(bar_to_slit(1) == 1)
        sa(bar_to_slit(2) == 1)
        sa(bar_to_slit(91)==46)
        sa(bar_to_slit(92)==46)
        sa(bar_to_slit(92.)==46)
        sa(bar_to_slit(1.)==1)
        sa(bar_to_slit(1.5)==1)
        sa(bar_to_slit(2.)==1)
        self.assertRaises(MismatchError, bar_to_slit, (-1, 0, 93, 94))

    def test_bar_mm(self):
        sa = self.assertTrue


        # Values are taken from ccs
        #p0 = mosfire_geoxytran(0,0)
        p0 = python_geoxytran(0,0)
        sa(np.abs(p0[0] - center_pix[0]) < 1e-6)
        sa(np.abs(p0[1] - center_pix[1]) < 1e-6)

    def test_Barset(self):
        b = Barset()
        pos = np.arange(92)
        b.set_mms(pos)
        self.assertTrue((b.pos == pos).all())

        p1 = b.get_bar_pix(1)
        #self.assertTrue((p1==csu_mm_to_pix(pos[0], 1)).all())

        # FIXME More tests here

    
if __name__ == '__main__':
    unittest.main()
