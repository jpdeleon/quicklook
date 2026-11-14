# Changelog

## [1.12.0](https://github.com/jpdeleon/quicklook/compare/v1.11.0...v1.12.0) (2026-11-14)


### Features

* give the single-event search its own wider detrending window ([111f089](https://github.com/jpdeleon/quicklook/commit/111f089c0ed43daebd43c19861d87717c5abd758))
* single-event search detrends with its own wider window ([#42](https://github.com/jpdeleon/quicklook/issues/42) step 3) ([933ad96](https://github.com/jpdeleon/quicklook/commit/933ad96238cebe7e7ecdfc83f17334390c29b0a1))
* TTV diagnostic byproduct figures for active stars (*_ttv.png, *_ttv_transits.png) ([b7dff94](https://github.com/jpdeleon/quicklook/commit/b7dff94bba186c6bcda801ced7581e57d884bf16))
* TTV diagnostic byproduct figures from local-baseline timing on raw flux ([25d5254](https://github.com/jpdeleon/quicklook/commit/25d52541c1a9b7ac933d7bd5ee2bd3eef4e4d008))
* TTV search with per-transit timing fits and O-C significance ([ace7e32](https://github.com/jpdeleon/quicklook/commit/ace7e320bc161ae5f904409f899eeb3be3aa2dfd))
* TTV search with per-transit timing fits and O-C significance ([#42](https://github.com/jpdeleon/quicklook/issues/42) step 5) ([5aec6cb](https://github.com/jpdeleon/quicklook/commit/5aec6cb6470b9fc823ff4249030f94c60b357c9e))


### Bug Fixes

* --custom-ephem crashed with "'tuple' object does not support item ([25d5254](https://github.com/jpdeleon/quicklook/commit/25d52541c1a9b7ac933d7bd5ee2bd3eef4e4d008))

## [1.11.0](https://github.com/jpdeleon/quicklook/compare/v1.10.1...v1.11.0) (2026-11-13)


### Features

* add opt-in single-transit (mono) search with none/mono/periodic classification ([809dbd0](https://github.com/jpdeleon/quicklook/commit/809dbd023f3cd792b503aaed25568016767c7d3f))
* break the time axis between separate sectors in mono/duo panels ([6b9646f](https://github.com/jpdeleon/quicklook/commit/6b9646f7a34f3e2a7478d42239d9e374914437d6))
* pair single events into duotransits with ranked period aliases ([6880468](https://github.com/jpdeleon/quicklook/commit/6880468cfd3bc4ed1f6c27e09e739be363ab0f62))
* pair single events into duotransits with ranked period aliases ([e8e06d0](https://github.com/jpdeleon/quicklook/commit/e8e06d0d05bb2fb3a51d8884e70e39bbb0e4ea3b))
* single-transit (mono) search with none/mono/periodic classification ([bfd7bd4](https://github.com/jpdeleon/quicklook/commit/bfd7bd44a1479d1e04ca8bdadf7eae4bb972c178))

## [1.10.1](https://github.com/jpdeleon/quicklook/compare/v1.10.0...v1.10.1) (2026-11-13)


### Performance Improvements

* cache TLS summary rows so /tls-summary stops hanging the GUI ([2d42a47](https://github.com/jpdeleon/quicklook/commit/2d42a47e32b0c573abcf6e69e30321965cfa1cf7))
* cache TLS summary rows so /tls-summary stops hanging the GUI ([ac45977](https://github.com/jpdeleon/quicklook/commit/ac4597755bdc9307acd3d3b011ef684adf2d8c16))

## [1.10.0](https://github.com/jpdeleon/quicklook/compare/v1.9.0...v1.10.0) (2026-11-12)


### Features

* add TOI-519 sectors 7 and 8 to the TGLC benchmark targets ([fcdff17](https://github.com/jpdeleon/quicklook/commit/fcdff175570c6b9b334165e9e167206e9d865a99))
* benchmark local TGLC light curves against the MAST HLSP ([57ba4f3](https://github.com/jpdeleon/quicklook/commit/57ba4f36c3d4f380aed4ed6c29e0282c22c30010))
* export TGLC quality flags and mask them in ql's local TGLC path ([6c783d7](https://github.com/jpdeleon/quicklook/commit/6c783d7a446938053b168ec415098e52d1d5cb00))
* expose the field-star prior in get_tglc_lc and the TGLC benchmark ([728f6f8](https://github.com/jpdeleon/quicklook/commit/728f6f803c17220740cb3e12c957e0684df2b4c3))


### Bug Fixes

* extract the requested target in get_tglc_lc, not the brightest star ([d1f24ac](https://github.com/jpdeleon/quicklook/commit/d1f24acd02c6ad7a7e171413ad3d63a7af93e7a1))

## [1.9.0](https://github.com/jpdeleon/quicklook/compare/v1.8.0...v1.9.0) (2026-11-11)


### Features

* extend TLS period ceiling to ~baseline and add catalogue epoch check ([242dd97](https://github.com/jpdeleon/quicklook/commit/242dd97e1424246aee03a804154163b53fb0ef36))
* extend TLS period ceiling to ~baseline and add catalogue epoch check ([c7f0db9](https://github.com/jpdeleon/quicklook/commit/c7f0db99de7a15e8f95b1953cec5bd37a7f10613))

## [1.8.0](https://github.com/jpdeleon/quicklook/compare/v1.7.0...v1.8.0) (2026-11-09)


### Features

* neighbor_check option; save Prot_ls and clear failed GLS ([1dc5b1f](https://github.com/jpdeleon/quicklook/commit/1dc5b1f76ccf296edd5f5053aa8707cb2a5c5f46))
* six detection/vetting improvements from gallery false-positive review ([7773344](https://github.com/jpdeleon/quicklook/commit/77733443dcdcd0e970ca416b7794f17bc5624355))


### Bug Fixes

* guard against missing per_transit_count in read-tls summary ([102b089](https://github.com/jpdeleon/quicklook/commit/102b08936b044a584b951cb894e2068a1ec9900a))
* penalize depth_ratio symmetrically in select_best_window ([933df68](https://github.com/jpdeleon/quicklook/commit/933df682350e426b0a8333b7b1cf69d25883e76a))

## [1.7.0](https://github.com/jpdeleon/quicklook/compare/v1.6.0...v1.7.0) (2026-11-05)


### Features

* add --sectors to stitch an explicit set of sectors ([ee3af86](https://github.com/jpdeleon/quicklook/commit/ee3af86038a0ba892a44245da1db14d317a2559b))
* add --sectors to stitch an explicit set of sectors ([508801a](https://github.com/jpdeleon/quicklook/commit/508801ab8025421daeb85e08ebf5719cb36c77ab))
* support --sector all in the run CLI ([dd853f1](https://github.com/jpdeleon/quicklook/commit/dd853f1fc81ad89636b30c5aa44615a0074d503b))
* support --sector all in the run CLI ([cb6ff3d](https://github.com/jpdeleon/quicklook/commit/cb6ff3df47da34e1e393552f203413faea52f2be))


### Bug Fixes

* do not abort run when Simbad query fails ([472b5f3](https://github.com/jpdeleon/quicklook/commit/472b5f378893ce603dc5055c9c936f7c85a12f07))
* do not abort run when Simbad query fails ([3948602](https://github.com/jpdeleon/quicklook/commit/39486023dfb5d201e35069603bcda45cf7655c49))
