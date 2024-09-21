# syntax=docker/dockerfile:1

# build gerris package
FROM debian:bookworm AS builder

ENV DEBIAN_FRONTEND=noninteractive

# enable sources
COPY <<EOF /etc/apt/sources.list.d/debian-src.sources
Types: deb-src
URIs: http://deb.debian.org/debian/
Suites: stable stable-updates
Components: main
Signed-By: /usr/share/keyrings/debian-archive-keyring.gpg
EOF

WORKDIR /build 

RUN <<EOF
apt-get update
apt-get build-dep -y gerris
apt-get source gerris
EOF

COPY <<EOF /build/gerris-nrrd.patch
diff -ur ./src/output.c /home/garth/gerris/src/output.c
--- ./src/output.c	2013-12-06 16:51:21.000000000 +0100
+++ /home/garth/gerris/src/output.c	2014-01-29 09:35:27.963214383 +0100
@@ -1425,6 +1425,18 @@
       break;
     }
 
+    case GFS_NRRD: {
+      gfs_domain_write_nrrd (domain, output->max_depth, domain->variables_io, output->precision,
+			    GFS_OUTPUT (event)->file->fp);
+      break;
+    }
+
+    case GFS_BOV: {
+      gfs_domain_write_bov (domain, output->max_depth, domain->variables_io, output->precision,
+			    GFS_OUTPUT (event)->file->fp);
+      break;
+    }
+
     case GFS_TECPLOT: {
       gfs_domain_write_tecplot (domain, output->max_depth, domain->variables_io, output->precision,
 				GFS_OUTPUT (event)->file->fp);
@@ -1474,6 +1486,7 @@
   switch (output->format) {
   case GFS_TEXT:    fputs (" format = text", fp);    break;
   case GFS_VTK:     fputs (" format = VTK", fp);     break;
+  case GFS_NRRD:    fputs (" format = NRRD", fp);    break;
   case GFS_TECPLOT: fputs (" format = Tecplot", fp); break;
   default: break;
   }
@@ -1539,6 +1552,8 @@
 	output->format = GFS_TEXT;
       else if (!strcmp (format, "VTK"))
 	output->format = GFS_VTK;
+      else if (!strcmp (format, "NRRD"))
+	output->format = GFS_NRRD;
       else if (!strcmp (format, "Tecplot"))
 	output->format = GFS_TECPLOT;
       else {
diff -ur ./src/output.h /home/garth/gerris/src/output.h
--- ./src/output.h	2013-12-06 16:51:21.000000000 +0100
+++ /home/garth/gerris/src/output.h	2014-01-29 09:35:27.964214398 +0100
@@ -147,6 +147,8 @@
 typedef enum   { GFS, 
 		 GFS_TEXT, 
 		 GFS_VTK, 
+		 GFS_NRRD,
+		 GFS_BOV,
 		 GFS_TECPLOT }              GfsOutputSimulationFormat;
 
 struct _GfsOutputSimulation {
diff -ur ./src/unstructured.c /home/garth/gerris/src/unstructured.c
--- ./src/unstructured.c	2013-12-06 16:51:21.000000000 +0100
+++ /home/garth/gerris/src/unstructured.c	2014-01-29 09:35:28.037215448 +0100
@@ -305,6 +305,241 @@
     gts_object_destroy (GTS_OBJECT (v[i]));
 }
 
+static void write_raw_cell(FttCell * cell, gpointer * param)
+{
+  FttVector* min    = param[1];
+  FttVector* max    = param[2];
+  gint*      isize  = param[3];
+  gfloat*    data   = param[4];
+  FttVector* lambda = param[5];
+  gdouble    csize  = ftt_cell_size(cell)/2.;
+
+  FttVector p, p1, p2;
+  gint i1, j1, i2, j2, k1, k2, i, j, k;
+
+  ftt_cell_pos (cell, &p);
+
+  p.x = p.x - csize / lambda->x + 1e-9;
+  p.y = p.y - csize / lambda->y + 1e-9;
+
+  p1.x = (p.x - csize)/lambda->x + 1e-9;
+  p1.y = (p.y - csize)/lambda->y + 1e-9;
+  p1.z = (p.z - csize)/lambda->z + 1e-9;
+  p2.x = (p.x + csize)/lambda->x - 1e-9;
+  p2.y = (p.y + csize)/lambda->y - 1e-9;
+  p2.z = (p.z + csize)/lambda->z - 1e-9;
+
+  // image_draw_square (image, &p1, &p2, c);
+  i1 = (p1.x - min->x)/(max->x - min->x)*isize[0];
+  i2 = (p2.x - min->x)/(max->x - min->x)*isize[0];
+  j1 = (p1.y - min->y)/(max->y - min->y)*isize[1];
+  j2 = (p2.y - min->y)/(max->y - min->y)*isize[1];
+  k1 = (p1.z - min->z)/(max->z - min->z)*isize[2];
+  k2 = (p2.z - min->z)/(max->z - min->z)*isize[2];
+
+  for( i=i1; i<=i2; ++i )
+  {
+    for( j=j1; j<=j2; ++j )
+    {
+      for( k=k1; k<=k2; ++k )
+      {
+          GfsVariable* v = param[6];
+          data[k*isize[0]*isize[1]+j*isize[0]+i] = GFS_VALUE(cell, v);
+      }
+    }
+  }
+}
+
+
+static void max_extent (FttCell * cell, FttVector * max)
+{
+  FttVector pos;
+
+  ftt_cell_pos (cell, &pos);
+  if (pos.x > max->x) max->x = pos.x;
+  if (pos.y > max->y) max->y = pos.y;
+  if (pos.z > max->z) max->z = pos.z;
+}
+
+
+static void min_extent (FttCell * cell, FttVector * min)
+{
+  FttVector pos;
+
+  ftt_cell_pos (cell, &pos);
+  if (pos.x < min->x) min->x = pos.x;
+  if (pos.y < min->y) min->y = pos.y;
+  if (pos.z < min->z) min->z = pos.z;
+}
+
+
+
+/**
+ * gfs_domain_write_nrrd:
+ * @domain: a #GfsDomain.
+ * @max_depth: the maximum depth to consider.
+ * @variables: a list of #GfsVariable to output.
+ * @precision: the formatting string for converting float to ASCII.
+ * @fp: a file pointer.
+ *
+ * Writes in @fp a NRRD-formatted representation of @domain and of the
+ * corresponding variables in the given list.
+ */
+void gfs_domain_write_nrrd (GfsDomain * domain, gint max_depth, GSList * variables,
+			   const gchar * precision, FILE * fp)
+{
+  gint size;
+  gint depth;
+  gint isize[3];
+  gdouble h;
+  gpointer param[8];
+  gfloat* varptr;
+  guint i;
+  guint *uiptr;
+  unsigned char* ucptr;
+  GSList* vi = variables;
+
+  FttVector cmax = { -G_MAXDOUBLE, -G_MAXDOUBLE, -G_MAXDOUBLE };
+  FttVector cmin = { G_MAXDOUBLE, G_MAXDOUBLE, G_MAXDOUBLE };
+
+  g_return_if_fail (domain != NULL);
+  g_return_if_fail (precision != NULL);
+  g_return_if_fail (fp != NULL);
+
+  if( max_depth < 0 )
+    size = 1 << gfs_domain_depth(domain);
+  else
+    size = 1 << max_depth;
+
+  gfs_domain_cell_traverse (domain, FTT_PRE_ORDER, FTT_TRAVERSE_LEVEL,
+			                      domain->rootlevel, (FttCellTraverseFunc) min_extent, &cmin);
+  gfs_domain_cell_traverse (domain, FTT_PRE_ORDER, FTT_TRAVERSE_LEVEL,
+			                      domain->rootlevel, (FttCellTraverseFunc) max_extent, &cmax);
+
+  if( cmin.x == G_MAXDOUBLE )
+      return;
+
+  h = ftt_level_size (domain->rootlevel)/2.;
+
+  fprintf( fp,
+  	   "NRRD0001\\n"
+  	   "# Gerris simulation version %s (%s)\\n"
+  	   "#   variables ",
+  	   GFS_VERSION,
+  	   GFS_BUILD_VERSION );
+
+  vi = variables;
+
+  while( vi )
+  {
+      GfsVariable * v = vi->data;
+      fprintf( fp, "%s ", v->name );
+      vi = vi->next;
+  }
+
+  fprintf( fp, "\\n"
+	   "type: float\\n"
+	   "encoding: raw\\n"
+	   "endian: little\\n" );
+
+#if FTT_2D
+  cmin.x = (cmin.x - h)/domain->lambda.x;
+  cmax.x = (cmax.x + h)/domain->lambda.x;
+  cmin.y = (cmin.y - h)/domain->lambda.y;
+  cmax.y = (cmax.y + h)/domain->lambda.y;
+
+  isize[0] = (cmax.x - cmin.x) * size;
+  isize[1] = (cmax.y - cmin.y) * size;
+  isize[2] = 1;
+
+  fprintf( fp,
+	   "dimension: 3\\n"
+	   "sizes: %d %d %d\\n"
+	   "spacings: %f %f NaN\\n"
+	   "centers: cell cell none\\n",
+	   isize[0],
+	   isize[1],
+	   g_slist_length(variables),
+	   (cmax.x - cmin.x)/isize[0],
+	   (cmax.y - cmin.y)/isize[1] );
+#else
+  cmin.x = (cmin.x - h)/domain->lambda.x;
+  cmax.x = (cmax.x + h)/domain->lambda.x;
+  cmin.y = (cmin.y - h)/domain->lambda.y;
+  cmax.y = (cmax.y + h)/domain->lambda.y;
+  cmin.z = (cmin.z - h)/domain->lambda.z;
+  cmax.z = (cmax.z + h)/domain->lambda.z;
+
+  isize[0] = (cmax.x - cmin.x) * size;
+  isize[1] = (cmax.y - cmin.y) * size;
+  isize[2] = (cmax.z - cmin.z) * size;
+
+  fprintf( fp,
+  	   "dimension: 4\\n"
+  	   "sizes: %d %d %d %d\\n"
+  	   "spacings: %f %f %f NaN\\n"
+  	   "centers: cell cell cell none\\n"
+	   "axis mins: %f %f %f NaN\\n"
+	   "axis maxs: %f %f %f NaN\\n",
+  	   isize[0],
+  	   isize[1],
+  	   isize[2],
+  	   g_slist_length(variables),
+  	   (cmax.x - cmin.x)/isize[0],
+  	   (cmax.y - cmin.y)/isize[1],
+  	   (cmax.z - cmin.z)/isize[2],
+	   cmin.x, cmin.y, cmin.z, 
+	   cmax.x, cmax.y, cmax.z );
+#endif
+
+  fprintf( fp, "\\n" );
+
+  varptr = malloc( sizeof(float)*isize[0]*isize[1]*isize[2] );
+
+  vi = variables;
+
+  while( vi )
+  {
+      GfsVariable* v = vi->data;
+
+      param[0] = fp;
+      param[1] = &cmin;
+      param[2] = &cmax;
+      param[3] = &isize;
+      param[4] = varptr;
+      param[5] = &domain->lambda;
+      param[6] = v;
+
+      memset( varptr, 0, sizeof(float)*isize[0]*isize[1]*isize[2] );
+
+      gfs_domain_cell_traverse( domain, FTT_PRE_ORDER, FTT_TRAVERSE_LEAFS|FTT_TRAVERSE_LEVEL, max_depth,
+                                (FttCellTraverseFunc) write_raw_cell, param );
+
+      fwrite( varptr, isize[0]*isize[1]*isize[2], sizeof(float), fp );
+
+      vi = vi->next;
+  }
+
+  free( varptr );
+}
+
+/**
+ * gfs_domain_write_bov:
+ * @domain: a #GfsDomain.
+ * @max_depth: the maximum depth to consider.
+ * @variables: a list of #GfsVariable to output.
+ * @precision: the formatting string for converting float to ASCII.
+ * @fp: a file pointer.
+ *
+ * Writes in @fp a BOV-formatted representation of @domain and of the
+ * corresponding variables in the given list.
+ */
+void gfs_domain_write_bov (GfsDomain * domain, gint max_depth, GSList * variables,
+			   const gchar * precision, FILE * fp)
+{
+  // TODO
+}
+
 static void write_tecplot_element (FttCell * cell, WriteParams * par)
 {
   static guint tecplot_index[NV] = {
diff -ur ./src/unstructured.h /home/garth/gerris/src/unstructured.h
--- ./src/unstructured.h	2013-12-06 16:51:21.000000000 +0100
+++ /home/garth/gerris/src/unstructured.h	2014-01-29 09:35:28.038215462 +0100
@@ -32,6 +32,16 @@
 			       GSList * variables, 
 			       const gchar * precision,
 			       FILE * fp);
+void gfs_domain_write_nrrd    (GfsDomain * domain,
+			       gint max_depth,
+			       GSList * variables,
+			       const gchar * precision,
+			       FILE * fp);
+void gfs_domain_write_bov     (GfsDomain * domain,
+			       gint max_depth,
+			       GSList * variables,
+			       const gchar * precision,
+			       FILE * fp);
 void gfs_domain_write_tecplot (GfsDomain * domain, 
 			       gint max_depth, 
 			       GSList * variables, 
EOF

WORKDIR /build/gerris-20131206+dfsg

RUN <<EOF
mv /build/gerris-nrrd.patch debian/patches
echo "gerris-nrrd.patch" >> debian/patches/series
dpkg-buildpackage -us -uc -j
EOF


FROM debian:bookworm

RUN --mount=type=bind,from=builder,source=/build,target=/build <<EOF 
apt-get update
apt install -y --no-install-recommends \
    /build/libgfs-1.3-2_20131206+dfsg-19_*.deb \
    /build/libgfs-dev_20131206+dfsg-19_*.deb \
    /build/gerris_20131206+dfsg-19_*.deb \
    m4 \
    python3-vtk9 \
    python3-pip \
    openmpi-bin
rm -rf /var/lib/apt/lists/*
EOF

RUN pip3 install  --break-system-packages pynrrd
