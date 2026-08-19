$(document).keyup(function(event) {
    if ($("#mz").is(":focus") && (event.keyCode == 13)) {
        $("#search").click();
    }
    if ($("#charge").is(":focus") && (event.keyCode == 13)) {
        $("#search").click();
    }
    if ($("#ppm").is(":focus") && (event.keyCode == 13)) {
        $("#search").click();
    }
});