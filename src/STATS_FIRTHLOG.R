#/***********************************************************************
# * Licensed Materials - Property of IBM 
# *
# * IBM SPSS Products: Statistics Common
# *
# * (C) Copyright IBM Corp. 2014
# *
# * US Government Users Restricted Rights - Use, duplication or disclosure
# * restricted by GSA ADP Schedule Contract with IBM Corp. 
# ************************************************************************/

# author__ = "SPSS, JKP"
# version__ = "1.0.2"

# History
# 22-apr-2014 Original Version
# 26-jun-2014 Adapt to incompatible changes in logistf vresion 1.20


helptext="STATS FIRTHLOG DEPENDENT=variable
    INDEP=list of variables
/OPTIONS FIRTH=YES or NO ALPHA=sig PPL=PROFILE or WALD
    MAXSTEP=number MAXHS=integer MAXIT=integer
    LCONV=number XCONV=number
/OUTPUT PLOT=list of integers
/SAVE DATASET=dataset name.

Compute Firth logistic regression.
Example:
STATS FIRTHLOG DEPENDENT=y INDEP=x1 x2 x3
/SAVE DATASET=casewiseresults.

DEPENDENT and INDEP list the dependent (binary) and
independent variables.  The dependent variable must
have a scale measurement level and have values 0 and 1.
All other parameters are optional.

FIRTH=YES specifies the use of Firth's penalized maximum likelihood
method.  NO specifies standard maximum likelihood.

PPL=PROFILE the use of the profile penalized log likelihood for
the confidence intervals and tests.  WALD specifies WALD tests.

CONF specifies the confidence level.  It must be a number between 
50 and 100.

MAXSTEP through XCONV specify iteration and convergence
criteria.  MAXIT is the maximum number of iterations, MAXSTEP
is the maximum step size in beta values per iteration.  MAXHS is the
maximum number of step halving in one iteration. LCONV is
the log likelihood convergence criterion and XCONV is the
change criterion for values in the beta vector.

PLOT specifies plots of the likelihood against the parameter value
for selected independent variables.  The variables are identified
by their number in the independent variables list starting from 1.
An index greater than the number of variables is ignored.  Categorical
variables are not plotted.

DATASET is the name for saving casewise results, which include
the predicted probability of the 1 value and the hat matrix
diagonal element for each case.  The first column is the
input case number.  The dataset is created even if
convergence is not achieved.

Notes:
Cases with missing data are always deleted.
SPSS weights are not honored.
With R3.1 and logistf 1.21, the hat matrix is not available.

STATS FIRTHLOG /HELP.  prints this information and does nothing else.
"

### MAIN ROUTINE ###
dofirth = function(dep, indep, firth=TRUE, conf=.95, ppl=TRUE, dataset=NULL,
        maxstep=NULL, maxhs=NULL, maxit=NULL, lconv=NULL, gconv=NULL, xconv=NULL,
        plotlist=NULL) {
    #estimate regression discontinuity model
    
    setuplocalization("STATS_FIRTHLOG")
    
    # A warnings proc name is associated with the regular output
    # (and the same omsid), because warnings/errors may appear in
    # a separate procedure block following the regular output
    procname=gtxt("Firth Logistic Regression")
    warningsprocname = gtxt("Firth Logistic Regression: Warnings")
    omsid="STATSFIRTH"
    warns = Warn(procname=warningsprocname,omsid=omsid)

    tryCatch(library(logistf), error=function(e){
        warns$warn(gtxtf("The R %s package is required but could not be loaded.", "logistf"),dostop=TRUE)
        }
    )
    logistfversion = as.numeric(packageDescription("logistf", fields="Version"))
    controllist = makecontrollist(list(maxit=maxit, maxhs=maxhs, maxstep=maxstep, lconv=lconv, 
                    gconv=gconv, xconv=xconv))
    alpha=(100. - conf)/100.
    frml = buildfrml(dep, indep)
    allargs = as.list(environment())
    if (!is.null(dataset)) {
        alldatasets = tolower(spssdata.GetDataSetList())
        if (tolower(dataset) %in% alldatasets) {
            warns$warn(gtxt("The output dataset name is already in use."), dostop=TRUE)
        }
        if ("*" %in% alldatasets) {
            warns$warn(gtxt("The active dataset must have a name to create an output dataset"), 
                dostop=TRUE)
        }
    }
    
    alldata = c(dep, indep)
    
    dta = spssdata.GetDataFromSPSS(alldata, missingValueToNA=TRUE, factorMode="levels")
    # validate the dependent variable
    if (is.factor(dta[[1]])) {
        warns$warn(gtxt("The dependent variable must have a scale measurement level"),
            dostop=TRUE)
    }
    depdist = table(dta[[1]], useNA="no")
    depvalues=attr(depdist, "dimnames")[[1]]
    if (length(depvalues) != 2 || length(intersect(depvalues, c(0,1))) != 2) {
        warns$warn(gtxt("The dependent variable has values other than 0 and 1 or is constant"),
            dostop=TRUE)
    }
    allargs["ndtarows"] = nrow(dta)
    # if saving a dataset, remove cases with missing values but record the case numbers
    if (!is.null(dataset)) {
        cc = complete.cases(dta)
        casenumbers = (1:nrow(dta))[cc]
        dta = dta[cc,]
    } else {
        casenumbers = NULL
    }
    res = tryCatch(logistf(formula=frml, data=dta, firth=firth, pl=ppl, 
            alpha=alpha, control=controllist),
        error = function(e) {
            warns$warn(e$message, dostop=TRUE)
            return(NULL)
        }
    )
    displayresults(allargs, dta, res, controllist, warns)
    if (!is.null(dataset)) {
        makedataset(allargs, res, casenumbers, warns)
    }
}

buildfrml = function(dep, indep) {
    # Return formula expression as formula object

    # dep is the name of dependent variable
    # indep is the list of names of the independent variables
    # warns is the error message object
    
    cov = paste(indep, collapse="+")

    frml = paste(dep, "~", cov, collapse=" ")
    return(as.formula(frml))
}

makecontrollist = function(alist) {
    # return control list with null items removed
    alist = alist[!sapply(alist, is.null)]
    return(do.call(logistf.control, alist))
}

displayresults = function(allargs, dta, res, controllist, warns) {
    # Produce pivot tables and charts
    
    StartProcedure(allargs[["procname"]], allargs[["omsid"]])
    summarylabels=list(
        gtxt("Dependent Variable"),
        gtxt("Conf. Interval Type"),
        gtxt("Conf. Interval Level (%)"),
        gtxt("Estimation Method"),
        gtxt("Output Dataset"),
        gtxt("Likelihood Ratio Test"),
        gtxt("Degrees of Freedom"),
        gtxt("Significance"),
        gtxt("Number of Complete Cases"),
        gtxt("Cases with Missing Data"),
        gtxt("Number of Iterations"),
        gtxt("Convergence Status"),
        gtxt("Last Log Likelihood Change"),
        gtxt("Maximum Last Beta Change")
    )

    lrt = -2*(res$loglik[1] - res$loglik[2])
    summaryvalues = list(
        allargs[["dep"]],
        ifelse(allargs[["ppl"]], gtxt("Profile penalized log likelihood"),
            gtxt("Wald")),
        allargs[["conf"]],
        ifelse(allargs[["firth"]], gtxt("Firth penalized maximum likelihood"),
            gtxt("maximum likelihood")),
        ifelse(is.null(allargs[["dataset"]]), gtxt("--NA--"), allargs[["dataset"]]),
        round(lrt, 4),
        res$df,
        1. - pchisq(lrt, res$df),
        res$n,
        allargs[["ndtarows"]] - res$n,
        res$iter,
        ifelse(res$conv[[1]] > controllist[["lconv"]] ||
            res$conv[[3]] > controllist[["xconv"]], gtxt("FAILED"), gtxt("Converged")),
        res$conv[[1]],
        res$conv[[3]]
    )

    names(summaryvalues) = summarylabels
    summarydf = data.frame(cbind(summaryvalues))
    colnames(summarydf) = gtxt("Values")

    spsspivottable.Display(summarydf, title=gtxt("Firth Logistic Regression Summary"), 
                           templateName="STATSFIRTHSUMMARY",
                           caption=gtxt("Results computed by R logistf package"),
                           isSplit=FALSE,
                           format=formatSpec.Count
    )
    
    ddf = data.frame(res$coef, sqrt(diag(res$var)), res$ci.lower, res$ci.upper,
        qchisq(1 - res$prob, 1), res$prob)
    names(ddf) = c(gtxt("Coefficient"), gtxt("Std. Error"), 
        gtxt("Lower CI"), gtxt("Upper CI"), gtxt("Chi-Square"), "Sig.")
    ddf[sapply(ddf, is.infinite)] = "."
    spsspivottable.Display(ddf, title=gtxt("Coefficients"),
        outline=gtxt("Coefficients"),
        templateName="STATSFIRTHCOEF", isSplit=FALSE,
        caption=gtxtf("Dependent Variable: %s", allargs[["dep"]])
    )

    if (!is.null(allargs[["plotlist"]])) {
        nindep = length(allargs[["indep"]])
        assign("allargs", allargs, envir=.GlobalEnv)
        # plots do not work for categorical predictors
        for (var in allargs[["plotlist"]]) {
            if (var <= nindep) {
                vname = allargs[["indep"]][[var]]
                if (!is.factor(dta[[vname]])) {
                    if (allargs$logistfversion < 1.20) {
                        tryCatch(
                        logistfplot(formula=allargs[["frml"]], data=dta,
                            which=as.formula(paste("~",vname, "-1")), firth=allargs[["firth"]],
                            alpha = allargs[["alpha"]], control=controllist),
                        error=function(e) {warns$warn(e$message, dostop=TRUE)}
                        )
                    }
                    else {
                        tryCatch(
                            plot(profile(res, variable=vname,
                            alpha=allargs[["alpha"]], control=controllist)),
                            error=function(e) {warns$warn(e$message, dostop=TRUE)}
                        )
                    }
                }
            }
        }
        rm("allargs", envir=.GlobalEnv)
    }
    spsspkg.EndProcedure()
}

makedataset = function(allargs, res, casenumbers, warns) {
    # create dataset containing predicted probabilities and hat values
    # casenumbers is used to account for missing data
    
    dict = list()
    dict[[1]] = c("Case", "", 0, "F8.0", "nominal")
    dict[[2]] = c("Prob", gtxt("Predicted probability"), 0, "F8.3", "scale")
    if (!is.null(res$hat.diag)) {
        dict[[3]] = c("Hat", gtxt("Hat Value"), 0, "F8.3", "scale")
    }
    dict = spssdictionary.CreateSPSSDictionary(dict)
    spssdictionary.SetDictionaryToSPSS(allargs[["dataset"]], dict)
    if (!is.null(res$hat.diag)) {
        spssdata.SetDataToSPSS(allargs[["dataset"]], data.frame(Case=casenumbers, 
            prob=res$predict, hat=res$hat.diag))
    } else {
        spssdata.SetDataToSPSS(allargs[["dataset"]], data.frame(Case=casenumbers, 
            prob=res$predict))  
    }
    spssdictionary.EndDataStep()
}

Warn = function(procname, omsid) {
    # constructor (sort of) for message management
    lcl = list(
        procname=procname,
        omsid=omsid,
        msglist = list(),  # accumulate messages
        msgnum = 0
    )
    # This line is the key to this approach
    lcl = list2env(lcl) # makes this list into an environment

    lcl$warn = function(msg=NULL, dostop=FALSE, inproc=FALSE) {
        # Accumulate messages and, if dostop or no message, display all
        # messages and end procedure state
        # If dostop, issue a stop.

        if (!is.null(msg)) { # accumulate message
            assign("msgnum", lcl$msgnum + 1, envir=lcl)
            # There seems to be no way to update an object, only replace it
            m = lcl$msglist
            m[[lcl$msgnum]] = msg
            assign("msglist", m, envir=lcl)
        } 

        if (is.null(msg) || dostop) {
            lcl$display(inproc)  # display messages and end procedure state
            if (dostop) {
                stop(gtxt("End of procedure"), call.=FALSE)  # may result in dangling error text
            }
        }
    }
    
    lcl$display = function(inproc=FALSE) {
        # display any accumulated messages as a warnings table or as prints
        # and end procedure state, if any

        if (lcl$msgnum == 0) {   # nothing to display
            if (inproc) {
                spss.EndProcedure()
            }
        } else {
            if (!inproc) {
                procok =tryCatch({
                    StartProcedure(lcl$procname, lcl$omsid)
                    TRUE
                    },
                    error = function(e) {
                        FALSE
                    }
                )
            }
            if (procok) {  # build and display a Warnings table if we can
                table = spss.BasePivotTable("Warnings ","Warnings") # do not translate this
                rowdim = BasePivotTable.Append(table,Dimension.Place.row, 
                    gtxt("Message Number"), hideName = FALSE,hideLabels = FALSE)

                for (i in 1:lcl$msgnum) {
                    rowcategory = spss.CellText.String(as.character(i))
                    BasePivotTable.SetCategories(table,rowdim,rowcategory)
                    BasePivotTable.SetCellValue(table,rowcategory, 
                        spss.CellText.String(lcl$msglist[[i]]))
                }
                spsspkg.EndProcedure()   # implies display
            } else { # can't produce a table
                for (i in 1:lcl$msgnum) {
                    print(lcl$msglist[[i]])
                }
            }
        }
    }
    return(lcl)
}

# localization initialization
setuplocalization = function(domain) {
    # find and bind translation file names
    # domain is the root name of the extension command .R file, e.g., "SPSSINC_BREUSCH_PAGAN"
    # This would be bound to root location/SPSSINC_BREUSCH_PAGAN/lang

    fpath = Find(file.exists, file.path(.libPaths(), paste(domain, ".R", sep="")))
    bindtextdomain(domain, file.path(dirname(fpath), domain, "lang"))
} 
# override for api to account for extra parameter in V19 and beyond
StartProcedure <- function(procname, omsid) {
    if (substr(spsspkg.GetSPSSVersion(),1, 2) >= 19) {
        spsspkg.StartProcedure(procname, omsid)
    }
    else {
        spsspkg.StartProcedure(omsid)
    }
}

gtxt <- function(...) {
    return(gettext(...,domain="STATS_FIRTHLOG"))
}

gtxtf <- function(...) {
    return(gettextf(...,domain="STATS_FIRTHLOG"))
}


Run = function(args) {
    #Execute the STATS FIRTHLOG command

    cmdname = args[[1]]
    args = args[[2]]
    oobj = spsspkg.Syntax(list(
        spsspkg.Template("DEPENDENT", subc="",  ktype="existingvarlist", var="dep"),
        spsspkg.Template("INDEP", subc="", ktype="existingvarlist", var="indep", islist=TRUE),
        
        spsspkg.Template("FIRTH", subc="OPTIONS", ktype="bool", var="firth"),
        spsspkg.Template("PPL", subc="OPTIONS", ktype="bool", var="ppl"),
        spsspkg.Template("PLOT", subc="OPTIONS", ktype="bool", var="doplot"),
        spsspkg.Template("CONF", subc="OPTIONS", ktype="float", var="conf",
            vallist=list(50,99.9999)),
        spsspkg.Template("MAXSTEP", subc="OPTIONS", ktype="float", var="maxstep"),
        spsspkg.Template("MAXHS", subc="OPTIONS", ktype="int", var="maxhs"),
        spsspkg.Template("MAXIT", subc="OPTIONS", ktype="int", var="maxit"),
        spsspkg.Template("LCONV", subc="OPTIONS", ktype="float", var="lconv"),
        spsspkg.Template("GCONV", subc="OPTIONS", ktype="float", var="gconv"),
        spsspkg.Template("XCONV", subc="OPTIONS", ktype="float", var="xconv"),
        
        spsspkg.Template("PLOT", subc="OUTPUT", ktype="int", var="plotlist", islist=TRUE,
            vallist=list(0)),
        
        spsspkg.Template("DATASET", sub="SAVE", ktype="varname", var="dataset")
    ))

    # A HELP subcommand overrides all else
    if ("HELP" %in% attr(args,"names")) {
        #writeLines(helptext)
        helper(cmdname)
    }
    else {
        res <- spsspkg.processcmd(oobj, args, "dofirth")
    }
}

helper = function(cmdname) {
    # find the html help file and display in the default browser
    # cmdname may have blanks that need to be converted to _ to match the file
    
    fn = gsub(" ", "_", cmdname, fixed=TRUE)
    thefile = Find(file.exists, file.path(.libPaths(), fn, "markdown.html"))
    if (is.null(thefile)) {
        print("Help file not found")
    } else {
        browseURL(paste("file://", thefile, sep=""))
    }
}
if (exists("spsspkg.helper")) {
assign("helper", spsspkg.helper)
}