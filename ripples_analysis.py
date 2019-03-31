import matplotlib
import matplotlib.pyplot as plt
import matplotlib.cm as cm
import matplotlib.gridspec as gridspec
import numpy as np
import cPickle as pickle
import glob


class Analysis(object):

    def __init__(self, xmlFile, burn=None):
        plt.close('all')
        self.xmlFile = xmlFile
        self.burn = burn
        self.read_par_file()

    def read_par_file(self):
        import xml.etree.ElementTree as ET
        tree = ET.parse(self.xmlFile)
        root = tree.getroot()

        ddict = dict([(child.tag, []) for child in root])

        pp = root[0]
        assert pp.tag == 'parameters'
        ppc = pp.getchildren()

        assert ppc[2].tag == 'zLENS'
        self.zLENS = ppc[2].text
        assert ppc[3].tag == 'zSOURCE'
        self.zSOURCE = ppc[3].text

        assert root[2].tag == 'action'
        aaa = root[2].getchildren()

        assert aaa[0].tag == 'Model'
        self.potential = aaa[0].text

        assert aaa[1].tag == 'parameter_flags'
        fixed = aaa[1].text
        self.fixed = list(map(int, fixed.split()))

        assert aaa[4].tag == 'task'
        self.task = aaa[4].text

        if self.task == 'MCMC':
            assert root[3].tag == 'MCMC'
            mmm = root[3].getchildren()
            assert mmm[0].tag == 'numWalkers'
            self.nwalker = mmm[0].text
            assert mmm[5].tag == 'output_filename'
            self.output_filename = mmm[5].text

        # read chain
        self.cFiles = glob.glob(self.output_filename + '/c*.txt')
        # print self.cFiles


    def get_chain(self, walker, column):
        cf = self.output_filename + '/chain_number_' + str(walker) + '.txt'
        chain = np.loadtxt(cf, usecols=(column))
        return chain


    def merge_chains_from_all_walkers(self, column):
        """

        Remove burn in steps before merging chain -- for plotting PDFs

        """

        for www in range(len(self.nwalker)):
            if www == 0:
                ccc = self.get_chain(www, column)
                if self.burn:
                    ccc = ccc[self.burn:]
            else:
                if self.burn:
                    _ccc = self.get_chain(www, column)[self.burn:]
                else:
                    _ccc = self.get_chain(www, column)
                ccc = np.hstack((ccc, _ccc))
        return ccc


    def plot_walker(self, saveFig=False):

        if self.potential == 'SIE':
            param_to_latex = dict(re=r"$R_E$",
                                  px=r"$P_x$",
                                  py=r"$P_y$",
                                  xcent=r"$x_{\rm cent}$",
                                  ycent=r"$y_{\rm cent}$")

            params = ["re", "px", "py", "xcent", "ycent"]
        else:
            raise NotImplementedError

        # plot
        # Create a figure object with same aspect ratio as a sheet of paper...
        fig = plt.figure(figsize=(16, 20.6))

        # I want the plot of individual walkers to span 2 columns
        gs = gridspec.GridSpec(len(params), 3)

        for ii, param in enumerate(params):

            ax1 = plt.subplot(gs[ii, :2])
            ax1.axvline(0,
                        color="#67A9CF",
                        alpha=0.7,
                        linewidth=2)

            for walker in np.arange(int(self.nwalker)):
                theChain = self.get_chain(walker, ii+1)  # first col is chi2
                ax1.plot(np.arange(len(theChain)) - self.burn, theChain,
                         drawstyle="steps",
                         alpha=0.5)

            ax1.set_ylabel(param_to_latex[param],
                           fontsize=16,
                           labelpad=18,
                           rotation="horizontal",
                           color='k')
            # Don't show ticks on the y-axis
            # ax1.yaxis.set_ticks([])

            # For the plot on the bottom, add an x-axis label. Hide all others
            if ii == len(params) - 1:
                ax1.set_xlabel("step number", fontsize=16,
                               labelpad=18, color='k')
            else:
                ax1.xaxis.set_visible(False)

            ax2 = plt.subplot(gs[ii, 2])

            # Create a histogram of all values past the burn steps.
            # Make 100 bins between the y-axis bounds defined by the 'walkers' plot.
            these_chains = self.merge_chains_from_all_walkers(ii+1)
            ax2.hist(np.ravel(these_chains),
                     bins=np.linspace(ax1.get_ylim()[0], ax1.get_ylim()[1], 100),
                     orientation='horizontal',
                     facecolor="#67A9CF",
                 edgecolor="none")

            # Same y-bounds as the walkers plot, so they line up
            ax1.set_ylim(np.min(these_chains), np.max(these_chains))
            ax2.set_ylim(ax1.get_ylim())
            ax2.xaxis.set_visible(False)
            ax2.yaxis.tick_right()

            # For the first plot, add titles and shift them up a bit
            if ii == 0:
                t = ax1.set_title("Walkers", fontsize=16, color='k')
                t.set_y(1.01)
                t = ax2.set_title("Posterior", fontsize=16, color='k')
                t.set_y(1.01)

            if param == "re" or param == "xcent" or param == "ycent":
                ax2.set_ylabel("arcsec",
                               fontsize=16,
                               rotation="horizontal",
                               color='k',
                               labelpad=18)
            ax2.yaxis.set_label_position("right")

            # Adjust axis ticks, e.g. make them appear on the outside of the plots and
            #   change the padding / color.
            ax1.tick_params(axis='x', pad=2, direction='out',
                            colors='k', labelsize=14)
            ax2.tick_params(axis='y', pad=2, direction='out',
                            colors='k', labelsize=14)

            # Removes the top tick marks
            ax1.get_xaxis().tick_bottom()

            # Hack because the tick labels for phi are wonky... but this removed the
            #   first and last tick labels so I can squash the plots right up against
            #   each other
            if param == "phi":
                ax2.set_yticks(ax2.get_yticks()[1:-2])
            else:
                ax2.set_yticks(ax2.get_yticks()[1:-1])

        if saveFig:
            fig.subplots_adjust(hspace=0.0, wspace=0.0, bottom=0.075,
                                top=0.9, left=0.12, right=0.88)
            plt.savefig(self.output_filename + '/walker.pdf', bbox_inches='tight')
        else:
            plt.show(block=False)


    def plot_chi2(self, saveFig=False):
        """
        cFiles = list of cFiles

        """

        plt.figure()
        for cff in self.cFiles:
            like = np.loadtxt(cff, usecols=0)
            plt.plot(like[self.burn:], alpha=0.9)

        plt.xlabel('steps')

        if saveFig:
            plt.savefig(self.output_filename + '/chi2.pdf', bbox_inches='tight')
        else:
            plt.show(block=False)


    def imaging_model(self, img, saveFig=True):

        plt.figure()
        plt.imshow(img, origin='lower')

        if saveFig:
            plt.savefig('lensed.pdf', bbox_inches='tight')
        else:
            plt.show(block=False)


    def plot_source(self, img, saveFig=True):

        plt.figure()
        plt.imshow(img, origin='lower')

        if saveFig:
            plt.savefig('src.pdf', bbox_inches='tight')
        else:
            plt.show(block=False)



if __name__ == '__main__':
    datapathPrefix='./data/RXJ0911/Analysis/'

    # from ripples_analysis import Analysis
    out = Analysis('./data/RXJ0911/pars_file_o.xml', burn=500)
    out.plot_chi2(saveFig=True)
    out.plot_walker(saveFig=True)

    srcIm = datapathPrefix + 'model__ModelSrc.txt'
    sss = np.loadtxt(srcIm)
    sss = sss.reshape((int(np.sqrt(len(sss))), -1))
    out.plot_source(sss, saveFig=True)

    # lensed
    lensedIm = np.loadtxt('MYim1.txt')
    lensedIm = lensedIm.reshape((int(np.sqrt(len(lensedIm))), -1))
    out.imaging_model(lensedIm, saveFig=True)
