#' Transform an adTree output into a QSM file, simplifies it and computes
#' tree volume and DBH
#'
#' @param mesh a mesh3d class object build from a .obj output from adTree.
#' @param min_diameter numeric. The minimum diameter to keep in the qsm and mesh
#'
#' @return a simplidfied mesh, a qsm object (i.e. a set of cylinders with
#'         startxyz, endxyz, ID, radius and length), the estimates for DBH, tree
#'         volume, tree total length and volume and legth by 5cm diameter class
#' @export
#'
#' @examples
#' \donttest{
#' # import the obj file output from adTree
#' file = system.file("extdata", "adtree_qsm.obj", package="lidUrb")
#' mesh = Morpho::obj2mesh(file)
#'
#' # plot the original mesh
#' rgl::open3d()
#' rgl::shade3d(mesh,col="black",add=T)
#'
#' # transform the adTree mesh into a qsm file, simplify the mesh and compute
#' # some woody structure features
#' qsm = lidUrb::adtree2qsm(mesh)
#'
#' # plot the simplified mesh
#' rgl::open3d()
#' rgl::shade3d(qsm$mesh,col="black",add=T)
#'
#' # plot simplified mesh with colored by cylinder ID
#' rgl::open3d()
#' rgl::shade3d(qsm$mesh,col=rep(qsm$qsm$cyl_ID,each=20),meshColor = "faces",add=T)
#'
#' # woodu structure features
#' qsm$DHB
#' qsm$Volume
#' qsm$Length
#' qsm$Diameter_class
#' }

adtree2qsm = function(mesh, min_diameter = 0.02){

  . = .GRP = ID = X = Y = Z = endX = endY = endZ = is.tip = it = meanX = meanY =
  meanZ = new_it = radius = radius_class = radius_cyl = startX = startY = startZ =
  volume = NULL

  if(class(mesh) != "mesh3d") stop("mesh must be of class mesh3d")
  if(!is.numeric(min_diameter)) stop("min_diameter must be numeric")

  # create an index for each cylinder in the mesh
  mesh_vec = rep(1:ncol(mesh$it),each = 20)

  # build a temporary table with the vertrices position and the cylinder ID
  temp = data.table::data.table(it = c(mesh$it), ID = rep(1:(ncol(mesh$it)/20),each=60))
  temp = unique(temp[order(temp$ID,temp$it)])
  # define if the vertices belong to the base or the tip of a cylinder
  temp[,is.tip := rep(c(rep(0,10),rep(1,10)),max(temp$ID))]

  # add the coordinates of each vertex
  temp[,':='(
    X = mesh$vb[1,it],
    Y = mesh$vb[2,it],
    Z = mesh$vb[3,it]
  )]

  # compute the mean starting and ending point (i.e. the tips of the segment
  # representing the center of the cylinder)
  temp[, ':='(
    meanX = mean(X),
    meanY = mean(Y),
    meanZ = mean(Z)
  ),by=.(ID,is.tip)]

  # compute the radius as the distance of the base center to the base points
  # end points not used due to varying shape of the tip
  temp[,radius := 0]
  temp[is.tip == 0,radius := sqrt( (X - meanX)^2 + (Y - meanY)^2 + (Z - meanZ)^2  )]
  temp[,radius := max(radius),by=ID]

  if(max(temp$radius) < min_diameter/2){
    stop("The diameter threshold is larger than the largest object in the tree. Please use a smaller threshold.")
  }

  # keep only the segments information
  temp = unique(temp[,.(meanX,meanY,meanZ,radius,ID)])

  # transform the temporary file into a QSM file that respects the aRchi format
  QSM = data.table::data.table(
    startX = temp$meanX[seq(1,(nrow(temp)-1),by = 2)],
    startY = temp$meanY[seq(1,(nrow(temp)-1),by = 2)],
    startZ = temp$meanZ[seq(1,(nrow(temp)-1),by = 2)],
    endX = temp$meanX[seq(2,nrow(temp),by = 2)],
    endY = temp$meanY[seq(2,nrow(temp),by = 2)],
    endZ = temp$meanZ[seq(2,nrow(temp),by = 2)],
    cyl_ID = temp$ID[seq(1,(nrow(temp)-1),by = 2)],
    radius_cyl = temp$radius[seq(1,(nrow(temp)-1),by = 2)]
  )
  rm(temp)

  # add length and volume for each cylinder
  QSM[,length := sqrt( (startX - endX)^2 + (startY - endY)^2 + (startZ - endZ)^2)]
  QSM[,volume := length*pi*radius_cyl^2]

  ############################################################
  # keep the part of the QSM with radius above given threshold
  ############################################################
  ##### in the QSM file
  QSM = QSM[radius_cyl >= min_diameter/2]

  ##### in the mesh
  # keep valid IDs
  mesh$it = mesh$it[,which(mesh_vec %in% unique(QSM$cyl_ID))]

  # define new vertices ID
  temp_it = data.table::data.table(it = c(mesh$it))
  temp_it[,new_it := .GRP,by = it]


  # keep only the parts of the mesh that belong the right radius class
  mesh$vb = mesh$vb[,unique(temp_it$it)]
  mesh$it = matrix(temp_it$new_it,nrow=3)
  mesh$normals = mesh$normals[,unique(temp_it$it)]

  rm(temp_it)

  # compute total tree volume
  total_volume = sum(QSM$volume)

  # compute total tree volume
  total_length = sum(QSM$length)

  # compute tree DBH
  DBH = QSM[ (startZ-min(startZ)) <= 1.3 & (endZ-min(startZ)) >= 1.3,radius_cyl]*2
  # if multistem : DBH is computed from the sum of the surface of the stems
  if(length(DBH > 1)){
    DBH = sqrt((sum(pi*(DBH/2)^2))/pi)*2
  }

  # compute the total length and volume by radius class (classes of 0.025 in radius, i.e. 5cm classes)
  QSM[,radius_class := ceiling(radius_cyl/0.025)*0.025]
  QSM[,radius_class := paste(format((radius_class*2)-0.05, nsmall = 2),format(radius_class*2, nsmall = 2),sep="-")]
  diam_class_metrics = QSM[,.(sum(length),sum(volume)),by = radius_class]
  data.table::setnames(diam_class_metrics,c("diameter_class","total_length","total_volume"))
  QSM[,radius_class := NULL]


  return(list(mesh = mesh,
              QSM = QSM,
              DBH = DBH,
              Volume = total_volume,
              Length = total_length,
              Diameter_class_metrics = diam_class_metrics))
}
