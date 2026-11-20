// Обращения к серверу. Ошибка сервера превращается в Error с текстом для пользователя.
async function unpack(resp) {
  let data = null;
  try { data = await resp.json(); } catch (_) { /* не JSON */ }
  if (!resp.ok) throw new Error((data && data.error) || `Сервер ответил с кодом ${resp.status}`);
  return data;
}

export function getJSON(path, params = {}) {
  const q = new URLSearchParams(params).toString();
  return fetch(q ? `${path}?${q}` : path).then(unpack);
}

export function postJSON(path, body) {
  return fetch(path, {
    method: "POST",
    headers: { "Content-Type": "application/json" },
    body: JSON.stringify(body),
  }).then(unpack);
}
