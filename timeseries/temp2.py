import geopandas as gpd

path=r"D:\PNM\PNM_crownmap_2025_fieldmaps_final\PNM_crownmap_2025_fieldmaps.shp"

gdf = gpd.read_file(path)

gdf.columns

gdf=gdf[['GlobalID', 'status',
       'latin', 'species_uk', 'notes', 'dead_stand', 'senecent_l', 'new',
       'Fruiting', 'Flowering', 'illuminati', 'liana', 'crown', 'tag',
       'leafing', 'EditDate', 'Editor', 'geometry']]


gdf= gdf.rename(columns={'GlobalID':'global_id',
                          'illuminati':'illumination',
                          'new': 'new_leaves',
                          'senecent_l': 'deciduous_leaves',
                          'Flowering': 'flowers',
                            'Fruiting': 'fruits',
                            'status': 'census_status',
                            'EditDate': 'census_date',
                            'species_uk': 'liana_species',
                            'Editor': 'collector',
                         })
gdf['area'] = gdf.area
gdf['plot']= "metrop" 

gdf['dead_stand'] = gdf['dead_stand'].fillna("no")
gdf['deciduous_leaves'] = gdf['deciduous_leaves'].fillna("no")
gdf['new_leaves'] = gdf['new_leaves'].fillna("no")
gdf['flowers'] = gdf['flowers'].fillna("no")
gdf['fruits'] = gdf['fruits'].fillna("no")
gdf['leafing'] = gdf['leafing'].fillna(100)

gdf.columns


gdf['flowers'] = gdf['flowers'].replace({'no': 'None', 'si': 'Partial'})
gdf['fruits'] = gdf['fruits'].replace({'no': 'None', 'si': 'Partial'})
gdf['new_leaves'] = gdf['new_leaves'].replace({'no': 'None', 'si': 'Partial'})
gdf['deciduous_leaves'] = gdf['deciduous_leaves'].replace({'no': 'None', 'si': 'Partial'})

gdf.to_file(r"D:\PNM\PNM_crownmap_2025.gpkg", driver='GPKG', index=False)
#need columns plot, area_m2, mnemonic, family, genus, specie, latin,


closeup= gpd.read_file(r"D:\PNM\PNM_2025_closeup_centroids_final\closeup_corrected.shp")
closeup.columns
closeup= closeup.rename(columns={'polygon_id':'image_id',
                                 'liana_pres':'liana_boolean',
                                    'liana_spp':'liana_species',
                                    'F_Label': 'filename',
                                    'GlobalID': 'global_id',
                                    'EditDate': 'census_date',
                                    'Editor': 'collector',
                                    'X_Longitud': 'X_Longitude'})
                                    
closeup=closeup[['image_id', 'status', 'liana_boolean', 'liana_species', 'notes',
       'latin', 'filename', 'X_Longitude', 'Y_Latitude', 'Z_Altitude',
       'global_id','census_date', 'collector',
       'geometry']]

closeup['liana_boolean'] = closeup['liana_boolean'].replace({'no': False, 'si': True})

closeup=closeup.to_crs(gdf.crs)
closeup['X_Longitude'] = closeup.geometry.x
closeup['Y_Latitude'] = closeup.geometry.y

#NOW this one is difficult

closeup['liana_species'] = closeup.apply(lambda row: gdf[gdf.geometry.intersects(row.geometry)]['liana_species'].iloc[0] if len(gdf[gdf.geometry.intersects(row.geometry)]) > 0 else None, axis=1)
closeup['global_id'] = closeup.apply(lambda row: gdf[gdf.geometry.intersects(row.geometry)]['global_id'].iloc[0] if len(gdf[gdf.geometry.intersects(row.geometry)]) > 0 else None, axis=1)
closeup['latin'] = closeup.apply(lambda row: gdf[gdf.geometry.intersects(row.geometry)]['latin'].iloc[0] if len(gdf[gdf.geometry.intersects(row.geometry)]) > 0 else None, axis=1)

closeup.to_file(r"D:\PNM\PNM_crownmap_2025_closeup.gpkg", driver='GPKG', index=False)