

#' Compute topology from a qsm file
#'
#' @param QSM_file a qsm file
#'
#' @return the qsm file now includin axid_ID, axis sections and branching order
#' @export
#'
#' @examples
#' \donttest{
#' # import the obj file output from adTree
#' file = system.file("extdata", "adtree_qsm.obj", package="lidUrb")
#' mesh = Morpho::obj2mesh(file)
#'
#' # transform the adTree mesh into a qsm file, simplify the mesh and compute
#' # some woody structure features
#' qsm = lidUrb::adtree2qsm(mesh)
#'
#' # compute qsm topology
#' qsm$qsm = lidUrb::qsm_topology(qsm$qsm)
#'
#' rgl::open3d()
#' rgl::shade3d(qsm$mesh,col=rep(qsm$qsm$branching_order,each=20),meshColor = "faces",add=T)
#'
#' rgl::open3d()
#' rgl::shade3d(qsm$mesh,col=rep(qsm$qsm$axis_ID,each=20),meshColor = "faces",add=T)
#' }
qsm_topology = function(QSM_file){

  . = .axis_ID = bear_length = branching_order = cyl_ID = parent_ID = section =
  axis_ID = volume = NULL

  ###################################
  ##### compute tree topology #######
  ###################################

  QSM = QSM_file$QSM

  QSM[,cyl_ID := 1:nrow(QSM)]

  # compute parent cylinder ID
  QSM[,parent_ID := 0]
  QSM[,parent_ID := FNN::knnx.index(data = QSM[,4:6],query = QSM[,1:3], algorithm = "kd_tree", k = 1)]
  QSM[cyl_ID == parent_ID,parent_ID := 0]

  # compute total length bared by each cylinder
  QSM[, bear_length := 0]

  for(s in rev(sort(QSM$cyl_ID))){
    childs = QSM[cyl_ID == s | parent_ID == s]
    QSM[cyl_ID == s, bear_length := length+sum(childs$bear_length)]
  }


  ## the axis follows the longest bear_length
  QSM[, axis_ID := 0]

  cur_seg = QSM[parent_ID == 0] # start with the trunk base
  cur_ID = 1 # curent axis ID
  cur_sec = 1 # curent section (for segment computation)

  QSM[parent_ID == 0, axis_ID := cur_ID]

  queue = c()
  while(min(QSM$axis_ID)==0){
    childs = QSM[cyl_ID == cur_seg$cyl_ID,section := cur_sec]
    childs = QSM[parent_ID == cur_seg$cyl_ID]
    if(nrow(childs) >= 1){
      # if only one child -> it's in the same axis
      if(nrow(childs) == 1){
        QSM[cyl_ID == childs$cyl_ID, axis_ID := cur_ID]
        cur_seg = childs
      }
      # if more than one child -> the one that bare longest structure belongs
      # to the same axis, the other ones goes to the queue
      if(nrow(childs) > 1){
        QSM[cyl_ID == childs$cyl_ID[which.max(childs$bear_length)], axis_ID := cur_ID]
        cur_seg = childs[which.max(childs$bear_length),]
        queue = c(queue , childs$cyl_ID[-which.max(childs$bear_length)])

        cur_sec = cur_sec + 1
      }
    }else{
      # if there is no child -> pick the first one in the queue and increment cur_ID
      cur_ID = cur_ID + 1 # increment cur_ID
      cur_seg = QSM[cyl_ID == queue[1]] # select curent segment
      QSM[cyl_ID == cur_seg$cyl_ID, axis_ID := cur_ID]

      queue = queue[-1]

      cur_sec = cur_sec + 1
    }
  }

  # ADD BRANCHING ORDER
  cur_BO = QSM[QSM$axis_ID == 1] # axes of branching order 1
  QSM[axis_ID == 1,branching_order := 1]
  BO = 2 # first branching order to detect
  while(nrow(cur_BO)>0){
    # find all child axes of the axes of the curent BO
    child_axes=QSM[parent_ID %in% cur_BO$cyl_ID &
                     !(axis_ID %in% unique(cur_BO$axis_ID)),c(axis_ID)]
    # add the new BO to the child axes
    QSM[axis_ID %in% child_axes,branching_order := BO]
    # select the child axes for the next round
    cur_BO = QSM[QSM$axis_ID %in% child_axes]
    BO = BO+1
  }

  Branching_order_metrics = QSM[,.(sum(length),sum(volume)),by = branching_order]
  data.table::setnames(Branching_order_metrics,c("branching_order","total_length","total_volume"))

  Axes_metrics = QSM[,.(sum(length),sum(volume)),by = axis_ID]
  data.table::setnames(Axes_metrics,c("axis_ID","total_length","total_volume"))

  QSM[,':='(section = NULL, bear_length = NULL)]

  return(list(mesh = QSM_file$mesh,
              QSM = QSM,
              DBH = QSM_file$DBH,
              Volume = QSM_file$Volume,
              Length = QSM_file$Length,
              Diameter_class_metrics = QSM_file$Diameter_class_metrics,
              Branching_order_metrics = Branching_order_metrics,
              Axes_metrics = Axes_metrics)
         )
}





