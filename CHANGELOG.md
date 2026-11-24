# Changelog

## [1.14.0](https://github.com/jpdeleon/quicklook/compare/v1.13.0...v1.14.0) (2026-11-24)


### Features

* --precession-mode to run the EB precession branch from CLI and GUI ([8a6ccf6](https://github.com/jpdeleon/quicklook/commit/8a6ccf69918217c4f76d0ed2f25c0bd730506120))
* eclipsing-binary branch of the precession search ([#63](https://github.com/jpdeleon/quicklook/issues/63)) ([943e9a7](https://github.com/jpdeleon/quicklook/commit/943e9a7ff32c39d801c11f8d32dc2f61f65adb3c))
* eclipsing-binary branch of the precession search ([#63](https://github.com/jpdeleon/quicklook/issues/63)) ([175eeba](https://github.com/jpdeleon/quicklook/commit/175eebad7b46a101eabd45e0fab3125d224725d4))
* **gui:** add Precession search option and db/dt summary column ([26158c0](https://github.com/jpdeleon/quicklook/commit/26158c09bec477f1d11680d41a20b484a60b759d))
* precession search via per-sector transit shape ([#63](https://github.com/jpdeleon/quicklook/issues/63)) ([d1d3b02](https://github.com/jpdeleon/quicklook/commit/d1d3b02b43522f4c512917dadc3d6bbaf2b00512))
* search for nodal precession via per-sector transit shape ([#63](https://github.com/jpdeleon/quicklook/issues/63)) ([083d316](https://github.com/jpdeleon/quicklook/commit/083d3163169dc5f67eb82770880171fa78795d06))


### Bug Fixes

* block grazing precession trends only when Rp/Rs absorbs them ([3d3ac40](https://github.com/jpdeleon/quicklook/commit/3d3ac40ffdbfd4e4c7d9e217c9651387c62f3b51))
* exit 3 when the TGLC benchmark disagrees with MAST, keep 1 for crashes ([2f17435](https://github.com/jpdeleon/quicklook/commit/2f17435a9e1440a7539544b1467e83136a78a1fd))
* exit 3 when the TGLC benchmark disagrees with MAST, keep 1 for crashes ([eed08ef](https://github.com/jpdeleon/quicklook/commit/eed08ef0131d107d9d48ed31fcb25a79abe24cb6))
* **gui:** bound sector selections and preserve all-sector toggles ([13082c7](https://github.com/jpdeleon/quicklook/commit/13082c78c54302a7bf3f299ae75b5069e965801c))
* **gui:** display companion figures written by the selected job ([a00c510](https://github.com/jpdeleon/quicklook/commit/a00c510f3cfc2c5266d601a3ab41b0bb36b41b28))
* **gui:** multi-sector searches, result figures, and Colab notebooks for each search ([539793d](https://github.com/jpdeleon/quicklook/commit/539793dc2f508e9a3f34b3f0068fde7472d89d9d))
* **gui:** run multi-sector searches and show their companion figures ([59c0fe3](https://github.com/jpdeleon/quicklook/commit/59c0fe32f99ca484de0699c78e181fd01b44c697))
* load the optional precession fitter on demand ([0cc525e](https://github.com/jpdeleon/quicklook/commit/0cc525e8158051d7b0e847197377efc8f0223dec))
* normalize every stitched sector to the same baseline ([6c7577a](https://github.com/jpdeleon/quicklook/commit/6c7577a9b4f204dcbf39c4117497536a2376b024))
* preserve eclipsing-binary mode for skipped searches ([3e1067a](https://github.com/jpdeleon/quicklook/commit/3e1067a20bdaa3ff73534781144dc0dbe995a0dd))
* preserve known stellar errors and handle zero scatter ([fb28aa5](https://github.com/jpdeleon/quicklook/commit/fb28aa5748db387c21e8128ca50ba7efd2bac7c9))
* propagate joint impact-parameter covariance into trend tests ([911e0c4](https://github.com/jpdeleon/quicklook/commit/911e0c4068910f0f3ad9f47c8e59cb646cc54ed7))
* reject unconverged eclipse fits and test signed dilution trends ([0b6f143](https://github.com/jpdeleon/quicklook/commit/0b6f1435bf3b9fa5f0b428d05da5c24785670255))
* stop a linear ephemeris drift from raising t0_offset ([6601a07](https://github.com/jpdeleon/quicklook/commit/6601a07480c1766e28a3c80e5896802210b642e3))
* stop a missing distinct_transit_count from forcing periodic ([6a75833](https://github.com/jpdeleon/quicklook/commit/6a75833e31986f153df328f18580c31488d24052))

## [1.13.0](https://github.com/jpdeleon/quicklook/compare/v1.12.0...v1.13.0) (2026-11-15)


### Features

* broken time axes for multi-sector data in every layout ([360c509](https://github.com/jpdeleon/quicklook/commit/360c509f197915d4860f474741769dc1917668db))
* broken time axes for multi-sector data in every layout ([6d6c367](https://github.com/jpdeleon/quicklook/commit/6d6c3672c06b8c63b5976dd85922dffc991b44e5))
* ql shared-epochs flags single events shared across a batch ([3e57d37](https://github.com/jpdeleon/quicklook/commit/3e57d37b5f4e4ff12b9a8fda2e4f7e77865c1bec))
* red-noise timing errors and TTV vetting ([8a73c53](https://github.com/jpdeleon/quicklook/commit/8a73c531c2a372e0f840050b3b12c72ab4a8b0fe))
* red-noise TTV errors, TTV vetting, and cross-target shared-epoch check ([78acb67](https://github.com/jpdeleon/quicklook/commit/78acb6781c552f309b21732ca80d8b3faa45f1a2))


### Bug Fixes

* drop transits with degenerate timing errors ([304d319](https://github.com/jpdeleon/quicklook/commit/304d31983eaab26595b4bfc2137dc31960412bee))
* drop transits with degenerate timing errors ([7d171b5](https://github.com/jpdeleon/quicklook/commit/7d171b577771bd04239be32fe16123adbdee0c58))
* forward all options in --each-sector runs; rebuild TGLC caches lacking target coords ([62304e3](https://github.com/jpdeleon/quicklook/commit/62304e338ade2f5f1706ff166d20231a73b24576))
* forward every ql run option to --each-sector and --each-pipeline runs ([4b17c4a](https://github.com/jpdeleon/quicklook/commit/4b17c4a1a51faa7e8997c57df712e68db90e8f27))
* let --custom-ephem drive the TTV timing instead of TLS ([ba4752e](https://github.com/jpdeleon/quicklook/commit/ba4752e9119d184f4a5b435a001315b7ca0faea8))
* let --custom-ephem drive the TTV timing instead of TLS ([b282ef3](https://github.com/jpdeleon/quicklook/commit/b282ef361f684b1903473604bdb117377f02cee6))
* make GLS-based detrending tuning actually run ([992e5ff](https://github.com/jpdeleon/quicklook/commit/992e5ff0d161fa47946d03f5136f865602a1d72b))
* make GLS-based detrending tuning actually run ([b87547e](https://github.com/jpdeleon/quicklook/commit/b87547e8f568d7bb88be957994d14112e87b4143))
* rebuild TGLC source caches that lack target coordinates ([ab36077](https://github.com/jpdeleon/quicklook/commit/ab36077faa17839dd88fb063f798be75a891c286))
* refuse to detrend zero-centred or negative flux instead of hanging ([2985dd3](https://github.com/jpdeleon/quicklook/commit/2985dd3da62caa81ef296b79bd54a9ac56123eb0))
* refuse to detrend zero-centred or negative flux instead of hanging ([c9e6c39](https://github.com/jpdeleon/quicklook/commit/c9e6c39a9ce7c326a8bfa10a1ff4570b13356237))
* sort stitched multi-sector light curves by time ([d82e626](https://github.com/jpdeleon/quicklook/commit/d82e626cedaebf47506871cdab68b2f521d6ba36))
* sort stitched multi-sector light curves by time ([08e2c39](https://github.com/jpdeleon/quicklook/commit/08e2c39f083e97586b43d2081d5bd9a7f1116555))

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
* degenerate TLS fits and edge-of-CCD local TGLC extractions ([b87dd93](https://github.com/jpdeleon/quicklook/commit/b87dd9385a14812364118356a1222f409304a1cd))
* handle scalar per-transit TLS outputs when counting valid transits ([045e2e0](https://github.com/jpdeleon/quicklook/commit/045e2e01bddf0901cdc447e45eb69d0a5fd84a63))
* odd-even plot when no cadences fall in the transit window ([a5ac8f6](https://github.com/jpdeleon/quicklook/commit/a5ac8f6f875c9ee63ad7718881903b90298f0824))
* skip singular cycle-count fits in TTV search ([2cfdc94](https://github.com/jpdeleon/quicklook/commit/2cfdc948d0832b30f06b622a3834c12ae959d0fe))
* summary panel when neighbor-check panels are skipped ([21475cf](https://github.com/jpdeleon/quicklook/commit/21475cfe338d04eddf278dd2e4227610a7e74b8b))
* TGLC hardest bitmask and batch-run robustness ([f7f7fae](https://github.com/jpdeleon/quicklook/commit/f7f7faecda9e1eaa4c9a0df6755252c1bb779978))
* TGLC hardest bitmask and batch-run robustness ([2ff0d59](https://github.com/jpdeleon/quicklook/commit/2ff0d591267dfc109408d0a0da793cc2dcff35b2))

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
