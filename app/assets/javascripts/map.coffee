do (map = (window.map = window.map || {}), $ = jQuery) ->
  _map = null
  _markers = []
  _loading = null
  _timeout = 500

  initialize = () ->
    mapOptions =
      center: city
      zoom: zoom
    _map = new google.maps.Map(document.getElementById('map-canvas'), mapOptions)

    return

  google.maps.event.addDomListener(window, 'load', initialize)

  $('#create-meal-dialog').dialog
    autoOpen: false
    height: 300
    width: 600
    modal: true
    buttons:
      'Search': ()->
        query = $('#text-search-query').val()
        latitude = Number($('#text-search-latitude').val()).toFixed(6)
        longitude = Number($('#text-search-longitude').val()).toFixed(6)
        radius = Number($('#text-search-radius').val()).toFixed()

        if justWhiteSpace(query) or isNaN(latitude) or isNaN(longitude) or isNaN(radius) or radius < 0 or radius > 50000
          alert "Invalid data provided!"
          return false

        _text_search_query = query
        _text_search_latitude = latitude
        _text_search_longitude = longitude
        _text_search_radius = radius

        searchPlacesText(query, new google.maps.LatLng(latitude, longitude), radius)
        $(this).dialog('close')
        true
      'Cancel': ()->
        $(this).dialog('close')
        false

  $('#filter-hotel-list').button
    icons:
      primary: 'ui-icon-search'

  $('.save').button
    icons:
      primary: 'ui-icon-disk'
    text: false

  $('.save.name-button').button
    disabled: true

  $('#find-coordinates').button
    icons:
      primary: 'ui-icon-gear'

  $('#use-coordinates').button
    icons:
      primary: 'ui-icon-gear'

  $('#start-geocode').button
    icons:
      primary: 'ui-icon-search'


  $('.geocoder-link').click ()->
    w = window.open(geocoder_path, '_blank')
    w.focus()
    false

  $('#create-meal-button').click ()->
    $('#create-meal-dialog').dialog('open')
    false

  $('#zag-manual-info .save.chains').click ()->
    chain = $('input.chain').val()
    saveZagManualInfoChain(chain)
    false

  $('#zag-manual-info .save.addresses').click ()->
    address = $('input.address').val()
    city = $('input.city').val()
    state = $('input.state').val()
    zip = $('input.zip').val()
    country = $('input.country').val()
    formattedAddress = $('input.formattedAddress').val()
    saveZagManualInfoAddress(address, city, state, zip, country, formattedAddress)
    false


  $('#current-gds').click ()->
    marker = findMarkerForSource(_provider)
    $('#latitude').val(marker.latitude)
    $('#longitude').val(marker.longitude)
    true

  $('#text-search-current-gds-coordinates').click ()->
    marker = findMarkerForSource(_provider)
    $('#text-search-latitude').val(marker.latitude)
    $('#text-search-longitude').val(marker.longitude)
    true

  $('#find-coordinates').click ()->
    geocoder = new google.maps.Geocoder()
    request =
      address: $('#find-coordinates-input').val()
    geocoder.geocode(request, geocodeCallback)
    true

  $('#start-geocode').click ()->
    $.ajax
      url: geocoder_path
      type: 'post'
      data:
        provider: _provider
        identifier: _identifier
        name: _name
        chain: _chain
        country: _country
        city: _city
        address: _address
        iata: _iata
        use_coordinates: _use_coordinates
        latitude: _latitude
        longitude: _longitude
        radius: _radius
        no_zag_tag: _no_zag_tag
        zero_coordinates: _zero_coordinates
      beforeSend: () ->
        startLoading()
      error: (jqXHR, textStatus, errorThrown) ->
        stopLoading()
        alert "Error: " + textStatus
        return
      success: (data, textStatus, jqXHR) ->
        stopLoading()
        alert "Geocoding request was submitted"
        return
    true

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