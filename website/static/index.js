function deleteNote(cityId){
    fetch('/delete-note', {
        method: 'POST',
        body: JSON.stringify({ cityId: cityId}),
    }).then((_res) => {
        window.location.href ="/";
    });
}

function checkWeather(cityId){ 
    fetch('/weather',{
        method: 'POST',
        body: JSON.stringify({ cityId: cityId }),
    }).then((_res) => {
        window.location.href="/weather";   
    });
}