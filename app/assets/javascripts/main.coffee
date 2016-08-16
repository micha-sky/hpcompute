do (main = (window.main = window.main || {}), $ = jQuery) ->
  _map = null
  _markers = []
  _loading = null
  _timeout = 500
  _river = null;

  initialize = () ->
    mapOptions =
      center: city
      zoom: zoom
    _map = new google.maps.Map(document.getElementById('map-canvas'), mapOptions)

    $('#river-select').val('styr')
#    $('#river-select').selectmenu('refresh')

    return

  google.maps.event.addDomListener(window, 'load', initialize)

  $('#river-select').click () ->
    river = $(this).val()
    showRiverPoints(river)
    return

  showRiverPoints = (river) ->
    $.ajax
      method: 'get'
      url: 'get_points'
      data:
        river: river
      success: (points) ->
        clearAllMarkers()
        for point in points
          latLng = new google.maps.LatLng(point.latitude, point.longitude)
          marker = new google.maps.Marker(
            position: latLng,
            map: _map,
            icon: 'http://www.vitromolecularlaboratories.com/wp-content/themes/mercury/images/slider-dot.png'
          )
          marker.addListener 'click', () ->
            infowindow = buildInfoWindow(point)
            infowindow.open(_map, this)
          _markers.push(marker)
          setMapCenterAndZoom()
          console.log(point.point)
      error: (data) ->
        console.log(data)

  getResults = (river, branch, point) ->
    $.ajax
      type: 'get'
      url: 'get_results'
      data:
        river: river
        branch: branch
        point: point
      success: (data) ->
        console.log(data)
      error: (error) ->
        console.log(error)
    return

  clearAllMarkers = () ->
    marker.setMap(null) for marker in _markers
    _markers = []
    false

  buildInfoWindow = (point) ->
    contentString = '<div><b>River: ' + point.river +
    '</b></div>' + '<div><b>Branch: ' + point.branch +
    '</b></div>' + '<div><b>Point: ' + point.point + '</div>'

    infoWindow = new google.maps.InfoWindow(content: contentString)


  setMapCenterAndZoom = () ->
    bounds = getBoundsForMarkers()
    northEast = bounds.getNorthEast()
    southWest = bounds.getSouthWest()
    oneMarkerOnly = northEast.lat() == southWest.lat() and northEast.lng() == southWest.lng()

    if (oneMarkerOnly)
      _map.setCenter(bounds.getCenter())
      _map.setZoom(16)
    else
      _map.fitBounds(bounds)

    return

  getBoundsForMarkers = () ->
    bounds = new google.maps.LatLngBounds()
    for marker in _markers
      bounds.extend(marker.position)
    bounds

  $('#run-model').click ()->
    river = $('#river-select').val()
    showRiverPoints(river)

  runModel = (river) ->
    $.ajax
      url: run_model_path
      type: 'get'
      data:
        river: river
      beforeSend: () ->
        startLoading()
      complete: () ->
        stopLoading()
      success: (data) ->
        console.log(data)
      error: (error) ->
        console.log(error)
    return


  startLoading = () ->
    _loading = setTimeout("$('.loading').show()", _timeout)

  stopLoading = () ->
    clearTimeout(_loading)
    $('.loading').hide()

  clearAllMarkers = () ->
    marker.setMap(null) for marker in _markers
    _markers = []
    false

  clearAllMarkersZIndex = () ->
    for marker in _markers
      marker.infoWindow.setZIndex(1)
      marker.setZIndex(2)
    return

  getBoundsForMarkers = () ->
    bounds = new google.maps.LatLngBounds()
    for marker in _markers
      bounds.extend(marker.position)
    bounds

  searchPlacesText = (query, latLng, radius) ->
    request =
      location: latLng
      radius: radius
      query: query
      types: ['lodging']
    service = new google.maps.places.PlacesService(_map)
    service.textSearch(request, searchPlacesCallback)
    false

  searchPlacesNearby = (name, latLng, radius) ->
    request =
      location: latLng
      radius: radius
      name: name
      types: ['lodging']
    service = new google.maps.places.PlacesService(_map)
    service.nearbySearch(request, searchPlacesCallback)
    false

  searchPlacesCallback = (results, status) ->
    if status == google.maps.places.PlacesServiceStatus.OK
      addCandidateMarkers(results)
    else
      alert("Search failed!")
    false

  geocodeCallback = (results, status) ->
    if status == google.maps.GeocoderStatus.OK
      result = results[0]
      $('#find-coordinates-input').val(result.formatted_address)
      $('#latitude').val(result.geometry.location.lat().toFixed(6))
      $('#longitude').val(result.geometry.location.lng().toFixed(6))
    else
      alert("Google geocoding failed!")
    false

  textSearchGeocodeCallback = (results, status) ->
    if status == google.maps.GeocoderStatus.OK
      result = results[0]
      $('#text-search-find-coordinates-input').val(result.formatted_address)
      $('#text-search-latitude').val(result.geometry.location.lat().toFixed(6))
      $('#text-search-longitude').val(result.geometry.location.lng().toFixed(6))
    else
      alert("Google geocoding failed!")
    false

  setMapCenterAndZoom = () ->
    bounds = getBoundsForMarkers()
    northEast = bounds.getNorthEast()
    southWest = bounds.getSouthWest()
    oneMarkerOnly = northEast.lat() == southWest.lat() and northEast.lng() == southWest.lng()

    if (oneMarkerOnly)
      _map.setCenter(bounds.getCenter())
      _map.setZoom(16)
    else
      _map.fitBounds(bounds)

    return


  return