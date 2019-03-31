"""

Last mod: Dec 13 2018

Some util functions for preparing files and input parameters for Ripples.

To run in python (not CASA)

"""


def pick_highSNR(weight, bestNsam=30000, plothist=False):
    """
    Because of memory issues, should try to average more and select vis w/ lowest noise to start.

    Once you have a model, later you can run for a the full half a million or so visibilities by going on like 10-20 nodes.

    Hmm.. roblem with this though is that we are only selecting the shortest baselines..

    """

    import numpy as np
    if weight.ndim == 2:    # npol, nvis
        weight = np.average(weight, axis=0)

    elif weight.ndim > 2:
        raise NotImplementedError("Not implemnted. ")

    print "weight: "
    print "16, 50, 84th percentile: ", np.percentile(weight, (16, 50, 84))
    print "Max, Min: ", np.max(weight), np.min(weight)

    _sigma = 1./np.sqrt(weight)
    print "sigma: "
    print "16, 50, 84th percentile: ", np.percentile(_sigma, (16, 50, 84))
    print "Max, Min: ", np.max(_sigma), np.min(_sigma)


    if plothist:
        import matplotlib.pyplot as plt
        plt.figure()
        plt.hist(weight, bins=50)
        plt.show(block=False)

        plt.figure()
        plt.hist(weight, normed=True, cumulative=True, label='CDF',
                histtype='step', alpha=0.8, color='purple', bins=500)
        plt.show(block=False)


        plt.figure()
        plt.hist(_sigma, bins=50)
        plt.show(block=False)

        plt.figure()
        plt.hist(_sigma, normed=True, cumulative=True, label='CDF',
                histtype='step', alpha=0.8, color='purple', bins=500)
        plt.show(block=False)

    fract_of_pts = float(bestNsam) / len(weight)
    print fract_of_pts

    # threshold = np.percentile(_sigma, fract_of_pts*100)
    # indx, = np.where(_sigma < threshold)
    # # print "sigma and weight corresponding to cut: ", threshold, 1./threshold**2
    # assert len(_sigma[_sigma < threshold]) == bestNsam

    # ===========================================================
    # maybe this chuck is just wrong.......?
    # # to avoid only selecting the shortest baseline, where the sigma is the lowest..
    if 5.*fract_of_pts >= 0.8:
        if 2.5*fract_of_pts <= 0.8:
            fract_of_pts = 2.5*fract_of_pts
        else:
            pass
    else:
        fract_of_pts = 5.*fract_of_pts

    looseThreshold = np.percentile(_sigma, fract_of_pts*100)
    _indx, = np.where(_sigma < looseThreshold)
    # re-select from this pool of indexes
    indx = np.random.choice(_indx, size=bestNsam, replace=False)
    # ============================================================

    assert len(indx) == bestNsam

    return indx


def check_binFiles(path='./'):
    import glob, os
    import numpy as np

    ff = glob.glob(os.path.join(path, '*bin'))
    print ff
    for i in ff:
        xx = np.fromfile(i)
        print i, len(xx), xx.max(), xx.min()

    return None


def calc_img_from_vis(u, v, weight, vis, freq, imsize, pixsize, plot=True):
    """
    create image centered on zero.

    Parameters
    ----------
    u: array
        in meters
    v: array
        in meters
    weight: array
        length should be either twice that of u, v, with real_1, image_1, real_2, imag_2, real_3 ...., or same as that of u, v
    vis: array
        visibilities, length should be twice that of u, v, with real_1, imag_1, real_2, imag_2, etc
    freq: array
        array of frequencies in Hz, same length as u, v
    imsize: int
        size of output image in pix
    pixsize: float
        in arcsec

    Returns
    -------
    None

    """

    import time
    import numpy as np
    import matplotlib.pyplot as plt

    # from cy_modules.cython_helper import img_from_vis
    # time_now = time.time()
    # _zz = img_from_vis(u, v, weight, vis, imsize, pixsize)
    # time_end = time.time()
    # print time_end - time_now

    # time_now = time.time()
    arcsec2rad = np.pi / 180. / 3600.
    xxx = np.linspace(-imsize/2 * pixsize * arcsec2rad,
                        imsize/2 * pixsize * arcsec2rad,
                        imsize)
    X, Y = np.meshgrid(xxx, xxx)
    # print X, Y

    # CASA convension pixel alignment
    X -= pixsize / 2. * arcsec2rad
    Y += pixsize / 2. * arcsec2rad
    # print X, Y

    zz = np.zeros(X.shape).astype(complex)

    real = vis[0::2]
    imag = vis[1::2]

    # convert u, v [m] to lambda
    c = 299792458.0 # in m/s

    # print "u, v [m]", u, v
    lambda_ = c/freq   # m
    u /= lambda_
    v /= lambda_
    # print "u, v in lambda", u, v

    if len(weight) == len(u):
        for ii in range(len(u)):
            zz += weight[ii] * (real[ii] + 1j * imag[ii]) * np.exp(2j * np.pi * (u[ii] * X + v[ii] * Y))
    else:
        # weight in format of real_1, image_1, real_2, imag_2, real_3 ....
        for ii in range(len(u)):
            zz += weight[::2][ii] * (real[ii] + 1j * imag[ii]) * np.exp(2j * np.pi * (u[ii] * X + v[ii] * Y))
    # time_end = time.time()
    # print time_end - time_now

    if plot:
        plt.close('all')
        plt.figure()
        plt.subplot(111)
        plt.imshow(zz.real, origin='lower')
        plt.title('Real Part')
        # plt.subplot(122)
        # plt.imshow(zz.imag,  origin='lower')
        # plt.title('Imaginary Part')
        # plt.imshow(_zz.real,  origin='lower')
        # plt.title('Cython')
        try:
            plt.savefig('/mnt/ceph/users/interlens/Ripples/data/123.pdf')
            import os
            os.system('evince /mnt/ceph/users/interlens/Ripples/data/123.pdf & ')
        except:
            plt.show(block=False)

        # plt.savefig('/mnt/home/daisyleung/123.pdf')
        # import os
        # os.system('evince /mnt/home/daisyleung/123.pdf & ')

    # normalized to Jy/beam
    print zz.real.max(), zz.real.min()
    print (zz.real / (weight).sum()).max(), (zz.real / (weight).sum()).min()

    return zz.real / (weight).sum()


def test_image_vis(binpath='../', imsize=800, pixsize_arcsec=0.01, Nsam=None):
    import numpy as np
    import os

    # binary
    uu = np.fromfile(os.path.join(binpath, 'u.bin'))
    vv = np.fromfile(os.path.join(binpath, 'v.bin'))
    weight = np.fromfile(os.path.join(binpath, 'sigma_squared_inv.bin'))
    freqs = np.fromfile(os.path.join(binpath, 'frequencies.bin'))
    dataVis = np.fromfile(os.path.join(binpath, 'vis_chan_0.bin'))

    if Nsam:
        for _ in range(5):       # make 5 images, each with Npointss
            # already zipped the weight
            _idx = pick_highSNR(weight[::2], bestNsam=Nsam, plothist=False)
            _vis_real = [dataVis[iii*2: 2*iii+1][0] for iii in _idx]
            _vis_imag = [dataVis[iii*2+1: 2*iii+2][0] for iii in _idx]
            _vis = np.array(zip(_vis_real, _vis_imag)).flatten()

            im_jybeam = calc_img_from_vis(uu[_idx], vv[_idx], weight[
                          _idx * 2], _vis, freqs[_idx], imsize, pixsize_arcsec)
            import pdb; pdb.set_trace()
    else:
        im_jybeam = calc_img_from_vis(
            uu, vv, weight, dataVis, freqs, imsize, pixsize_arcsec)
    return im_jybeam


def plot_uv(u, v):
    import matplotlib.pyplot as plt
    plt.plot(u, v)
    plt.show(block=False)
    return None


def coord2pix(p_list, header):
    """
    return corresponding image pixel, given coordinate

    Utility:

    - to get lens center in pixel of image (corresponding to the input vis), based on coordinates
    - to get source pos in pixel of image


    Parameters
    ----------
    p_list: list of coords in astrolib.coords object

    Returns
    -------
    ra_px, dec_px

    """

    import pywcs
    ra_px_list = []
    dec_px_list = []
    wcs = pywcs.WCS(header).sub([1, 2])

    for p1 in p_list:
        ra, dec = p1
        ra_px, dec_px = wcs.wcs_sky2pix(ra, dec, 1)
        ra_px_list.append(ra_px[0])
        dec_px_list.append(dec_px[0])

    return ra_px_list, dec_px_list


def offset_from_phs(header, ra_px, dec_px):
    """
    Given imsize and ra and dec of a point in image pixel, return the x, y offset from the center of the image in pixel.

    If the spatial offset is large, we should do shiftpcd before extracting the visibilities as input for Ripples.

    """

    imsize = header['NAXIS1']
    center = int(imsize/2)

    x_cen, y_cen = ra_px - center, dec_px - center

    return x_cen, y_cen


def get_PB_arcsec(header):
    """ get PB size in arcsec"""

    import astropy.constants as c
    import numpy as np

    freqs = header['RESTFRQ']
    tele = header['TELESCOP']

    if 'ALMA' in tele:
        print("Assuming this is not ACA or TP data...")
        diam = 12.0   # m
    else:
        raise NotImplementedError("Not sure what's the dish size")

    clight = c.c.value
    PBfwhm = 1.2*(clight/freqs)/diam * (3600.*180./np.pi)
    print("PB FWHM: {:.2f} arcsec").format(PBfwhm)

    return None


def get_source_name(msFile):
    tb.open(msFile+'/SOURCE')
    objname = tb.getcol("NAME")
    objname = objname[0]
    tb.done()
    return objname



def get_nchan(lineFits):
    import pyfits
    import numpy as np

    fits_cube = pyfits.open(lineFits)
    header = fits_cube[0].header
    nchan = header['NAXIS3']

    return nchan, header



def vel2chan(lineFits, vel_kms, lineRestFreq_ghz=576.267931, redshift=2.7961):
    """
    Given velocity in km/s, return channel number in unbinned data..

    NOTE: potentially have issues if the line is broader than just one spw (with distinct frequencies -- okay if it's just multiple SPW of the same frequency tuning (e.g., multiple tracks))

    Parameters
    ----------
    lineFits: str
        filename of FITS

    vel_kms: float (or list of floats)
        velocity in km/s to get channel number

    Returns
    -------
    chan: channel number in the unbinned data corresponding to the input velocity


    """

    # import astropy.io.fits as fits
    import pyfits
    import numpy as np
    try:
        import pywcs
    except ImportError:
        print("Probably due to incompatibility in numpy version of CASA and pywcs. Re-run script in pure python instaed.")
    c_light = 299792.458    # km/s

    nchan, header = get_nchan(lineFits)
    chan_list = np.arange(nchan)

    vel_kms = list(vel_kms)
    wcs = pywcs.WCS(header)
    wcs_vel = wcs.sub([3])

    print header['CTYPE3'], header['CUNIT3']
    refFreq = header['CRVAL3']
    delFreq = header['CDELT3']

    freq_list = [refFreq + (i - header['CRPIX3']) * delFreq for i in chan_list ]
    freq = np.array(freq_list)

    restfreq = lineRestFreq_ghz / (1 + redshift)
    vel = (restfreq - freq/1.e9) / (restfreq) * c_light     # km/s

    # find nearest from vel
    chan = []
    for vvv in vel_kms:
        chan.append(np.argmin(np.abs(vel - vvv)))
    return chan


def shift_phs_center(vis, pcd, field):
    """
    Parameters
    ----------
    vis: str
        MS file
    pcd: str
        in format like: '10h27m51.6s -43d54m18s'
    field: str
    Returns
    -------
    outvis: str
        filename of the output vis, with the shifted phase center
    Example
    -------
    fixvis(vis=vis,
           outpuvis='ngc3256-fixed.ms',
           field='NGC3256',
           phasecenter='J2000 10h27m51.6s -43d54m18s')
    """

    outvis = vis.replace('.ms', 'phsShift.ms')
    fixvis(vis=vis, outputvis=outvis, field=field,
           phasecenter='J2000 ' + pcd)
    return outvis


def check_vis_fixvis(msFile, msFileFixvis):
    """
    compare visibilities from msFile before and after running fixvis()
    According to documentation, fixbis will recalculate the u,v,w coordinates relative to the new phase center.
    https://casa.nrao.edu/casadocs/casa-5.1.1/uv-manipulation/recalculation-of-uvw-values-fixvis
    """

    tb.open(msFile)
    data1 = tb.getcol('DATA')
    uvw1 = tb.getcol('UVW')
    tb.close()

    tb.open(msFileFixvis)
    data0 = tb.getcol('DATA')
    uvw0 = tb.getcol('UVW')
    tb.close()

    print data1[:100]
    print data0[:100]

    print uvw1[0, :100]
    print uvw0[0, :100]

    print uvw1[1, :100]
    print uvw0[1, :100]

    return None


def get_phs_center(msfile):
    """ get the phase center of the MS file using CASA taskinit.ms
    Parameters
    ----------
    msfile: str
        measurement set filename
    Returns
    -------
    pcd:
        phase center
    Note
    ----
    may need further modification if nField > 1 unless it is the same for all fields (need further investigation)
    """

    from taskinit import ms
    ms.open(msfile)
    pc = ms.getfielddirmeas()
    if not isinstance(pc, dict) is True:
        pc()
    epoch = pc['refer']
    pcd_dec = pc['m1']['value'] * 180 / np.pi
    pcd_ra = pc['m0']['value'] * 180 / np.pi
    if pcd_ra < 0:
        pcd_ra += 360
    ms.done()
    pcd = [pcd_ra, pcd_dec]
    return pcd



if __name__ == '__main__':

    import astrolib.coords as coords
    import pyfits

    # lens positions for RXJ0911, presumably we want to this to be the phase center of the data/input visibilities for Ripples!!!!
    p_list = [
         coords.Position("09:11:27.56447 05:50:54.302").j2000(),
         coords.Position("09:11:27.51420 05:50:54.967").j2000()
         ]

    fits_cube = pyfits.open('calibrated.LSRK_cont_R05.image.fits')
    header = fits_cube[0].header
    ra_px, dec_px = coord2pix(p_list, header)
    print "RA, Dec: "
    print zip(ra_px, dec_px)

    x_cen, y_cen = offset_from_phs(header, ra_px[0], dec_px[0])
    print("X, Y offset in px from ideal phase center: {:.2f}, {:.2f}").format(x_cen, y_cen)

    # for Ripples input
    get_PB_arcsec(header)

    # get channel number for slicing line data
    lineCubeFITS_unbinned = 'pseudo_line.contsub.nobin.R05.dirty.fits'
    vel_kms = [-100., 0.0]      # for test
    chan = vel2chan(lineCubeFITS_unbinned, vel_kms)
    print chan
    # update ripples_prep.yaml

