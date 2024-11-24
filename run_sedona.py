# run sedona on Polaris, ALCF
# 10/08/2024

import os
import geopandas as gpd
from pyspark.sql import SparkSession
from sedona.spark import *
import sys
import time


filename1 = sys.argv[1]
filename2 = sys.argv[2]


config = SedonaContext.builder() .\
    config('spark.jars.packages',
           'org.apache.sedona:sedona-spark-3.4_2.12:1.6.0,'
           'org.datasyslab:geotools-wrapper:1.6.0-28.2,'
           'uk.co.gresearch.spark:spark-extension_2.12:2.11.0-3.4'). \
            config('spark.executor.memory', '200g'). \
            config('spark.driver.memory', '200g'). \
    getOrCreate()

sedona = SedonaContext.create(config)


# Define rectangles as polygons using WKT format
# rectangles_data = [
#     (1, 'POLYGON((0 0, 0 2, 2 2, 2 0, 0 0))'),
#     (2, 'POLYGON((1 1, 1 3, 3 3, 3 1, 1 1))'),
#     (3, 'POLYGON((3.5 3, 3.5 5, 5 5, 5 3, 3.5 3))'),
#     (4, 'POLYGON((0.5 0.5, -0.5 0.5, -0.5 -0.5, 0.5 -0.5, 0.5 0.5))')
# ]

tick1 = time.time()

rectangle_data1 = []
with open(filename1) as file:
    for line in file:
        lst = line.strip().split()
        assert int(lst[0]) == 1
        rid = int(lst[1])
        low0 = float(lst[2])
        low1 = float(lst[3])
        high0 = float(lst[4])
        high1 = float(lst[5])
        rectangle_data1.append(
            (rid, 'POLYGON(({} {}, {} {}, {} {}, {} {}, {} {}))'.format(
                low0, low1,     \
                low0, high1,    \
                high0, high1,   \
                high0, low1,    \
                low0, low1
            ))
        )

tick2 = time.time()
print("loading first file takes: {}".format(tick2-tick1))


rectangle_data2 = []
with open(filename2) as file:
    for line in file:
        lst = line.strip().split()
        assert int(lst[0]) == 1
        rid = int(lst[1])
        low0 = float(lst[2])
        low1 = float(lst[3])
        high0 = float(lst[4])
        high1 = float(lst[5])
        rectangle_data2.append(
            (rid, 'POLYGON(({} {}, {} {}, {} {}, {} {}, {} {}))'.format(
                low0, low1,     \
                low0, high1,    \
                high0, high1,   \
                high0, low1,    \
                low0, low1
            ))
        )

tick3 = time.time()
print("loading second file takes: {}".format(tick3-tick2))


# Create Spark DataFrame with rectangles
rectangles_df1 = sedona.createDataFrame(rectangle_data1, ['id', 'geom'])
rectangles_df2 = sedona.createDataFrame(rectangle_data2, ['id', 'geom'])

# Register the DataFrame as a SQL temporary view
rectangles_df1.createOrReplaceTempView("rectangles1")
rectangles_df2.createOrReplaceTempView("rectangles2")

#sedona.sql("select * from rectangles1").show(5)
#sedona.sql("select * from rectangles2").show(5)

rectangles_with_geom1 = sedona.sql("""
SELECT id, ST_GeomFromText(geom) AS geom
FROM rectangles1
""")

rectangles_with_geom2 = sedona.sql("""
SELECT id, ST_GeomFromText(geom) AS geom
FROM rectangles2
""")

rectangles_with_geom1.printSchema()
rectangles_with_geom2.printSchema()

# Register the view with proper geometry data
rectangles_with_geom1.createOrReplaceTempView("rectangles_with_geom1")
rectangles_with_geom2.createOrReplaceTempView("rectangles_with_geom2")

# Query to find intersecting rectangles using ST_Intersects
intersecting_rectangles = sedona.sql("""
SELECT COUNT(a.id) AS intersections
FROM rectangles_with_geom1 a, rectangles_with_geom2 b
WHERE ST_Intersects(a.geom, b.geom)
""")

# Show the intersecting rectangles
# intersecting_rectangles.show()
print(intersecting_rectangles.show())


tick4 = time.time()
print("Find intersection takes: {}".format(tick4-tick3))
print("End to end takes: {}".format(tick4-tick1))